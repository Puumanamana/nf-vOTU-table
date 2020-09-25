nextflow.enable.dsl = 2

// Strategy
// 1) Merge contigs
// 2) Align reads against contigs
// 3) Normalize abundance table


// extract contigs matching the ids lists from multiple tools
process CONTIGS_FROM_IDS {
	tag {"${meta.id}"}
    container 'nakor/virus_extraction'
	publishDir "${params.outdir}/viral-fasta", mode: "copy"

	input:
    tuple val(meta), path(fasta), path(ids)

	output:
    tuple val(meta), path("viral_contigs-${meta.id}.fasta")
    
	script:
    def ids_all_tools = ids.join(' ')
    """
    contigs_from_ids.py \
        --fasta $fasta \
        --ids $ids_all_tools \
        --output viral_contigs-${meta.id}.fasta \
        --prefix "${meta.id}!"
    """
}

// concatenate and clusters the assemblies from multiple samples
// need to add a prefix (e.g. sample name) in case the contig ids
// are the same in the different assemblies
process MERGE_CONTIGS {
    publishDir "${params.outdir}/cluster-genomes", mode: 'copy'
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.fna'
    container 'nakor/votu_table'

    input:
    path f

    output:
    path '*', emit: all
    path '*.fna', emit: fasta
    path '*.clstr', emit: clusters    

    script:    
    """
    #!/usr/bin/env bash
    cat *.fasta > all_contigs.fasta
    Cluster_genomes.pl -f all_contigs.fasta -c 80 -i 95
    rm -f all_contigs.fasta
    """
}

process BWA_INDEX {
    container "nakor/votu_table"
    
    input:
    path assembly

    output:
    path 'ref*'

    script:
    """
    bwa index ${assembly} -p ref
    """
}

process BWA_MEM {
    container "nakor/votu_table"
    tag "${meta}" 
    publishDir "${params.outdir}/bwa", mode: 'copy', pattern: "*.{bam,bai}"
    
    input:
    tuple val(meta), path(fastq)
    path index

    output:
    path '*.bam', emit: bam
    path '*.bai', emit: bai
    
    script:
    """
    # Read alignment and filtering
    bwa mem -a -M -t $task.cpus ref $fastq \
        | samtools view -@ $task.cpus -bh -q ${params.bwa_qual} \
        | samtools view -@ $task.cpus -bh -f ${params.bwa_flag_f} -F ${params.bwa_flag_F} \
        | samtools sort -@ $task.cpus -o ${meta}.bam

    # Bam indexing
    samtools index -@ $task.cpus ${meta}.bam ${meta}.bai
    """
}


process MAKE_COVERAGE_TABLE {
    publishDir "${params.outdir}/bwa", mode: 'copy'
    container 'nakor/votu_table'

    input:
    path bam
    path bai

    output:
    path "coverage_table.tsv"

    script:
    def sorted_bam = bam.sort{it.getSimpleName()}
    def samples = sorted_bam.collect{it.getSimpleName()}
    """
    #!/usr/bin/env bash
    echo "contig pos ${samples.join(" ")}" | sed 's/ /\\t/g' > coverage_table.tsv
    samtools depth $sorted_bam >> coverage_table.tsv
    
    """
}

process MAKE_ABUNDANCE_TABLE {
    publishDir "${params.outdir}/tables", mode: 'copy'
    container 'nakor/votu_table'

    input:
    path table
    path fasta

    output:
    path "abundance_table.csv"

    script:
    """
    compute_table.py \
        --fasta ${fasta} \
        --coverage ${table}
    """
}


process SUBSET_ABUNDANCE_TABLE {
    publishDir "${params.outdir}/tables", mode: 'copy'
    container 'nakor/votu_table'

    input:
    val tools
    path viral_id_files
    path abundance_table
    path clusters

    output:
    path "abundance_table*.csv"

    script:
    """
    for tool in ${tools.join(' ')}; do
        subset_abundance_table.py \
            --name \$tool \
            --ids *\${tool}*.txt \
            --clusters ${clusters} \
            --abundance ${abundance_table}
    done
    """
}


workflow v_otu_table{
    take:
    contigs
    vir_ids
    reads

    main:
    // subset viral contigs
    // concatenates contigs from all samples
    viral_assembly = CONTIGS_FROM_IDS(
        contigs.combine(vir_ids, by: 0)
    )

    if (params.coassembly) {
        merged_assembly = viral_assembly.map{it[1]}
        clusters = file('empty')
    } else {        
        // merge assemblies from multiple samples
        merging = MERGE_CONTIGS(viral_assembly.collect{it[1]})
        merged_assembly = merging.fasta
        clusters = merging.clusters   
    }

    // align reads against clustered assembly
    ref = BWA_INDEX(merged_assembly).first() // to convert to value channel
    alignments = BWA_MEM(reads, ref)

    // compute coverage
    coverage_table = MAKE_COVERAGE_TABLE(
        alignments.bam.collect(),
        alignments.bai.collect()
    )

    // into abundance table
    abundance_table = MAKE_ABUNDANCE_TABLE(
        coverage_table,
        merged_assembly
    )

    tools = vir_ids
        .map{it[1]}
        .flatten()
        .map{"${it.getSimpleName().tokenize('-')[0]}"}
        .unique()
        .collect()

    SUBSET_ABUNDANCE_TABLE(
        tools,
        vir_ids.collect{it[1]},
        abundance_table,
        clusters
    )
}

workflow {
    // contig id are sample names
    contigs = Channel.fromPath(params.fasta).map{[[id: "${it.getSimpleName().tokenize('-')[-1]}"], it]}
    // vir_ids filename are {tool}-XXX-{sample_name}.txt
    vir_ids = Channel.fromPath(params.ids)
        .map{[[id: "${it.getSimpleName().tokenize('-')[-1]}"], it]}
        .groupTuple() // emits tuple(sample, files)
    reads = Channel.fromFilePairs(params.fastq)

    v_otu_table(contigs, vir_ids, reads)
}
