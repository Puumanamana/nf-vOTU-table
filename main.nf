nextflow.enable.dsl = 2

// Strategy
// 1) Merge contigs
// 2) Align reads against contigs
// 3) Normalize abundance table


// extract contigs matching the ids lists from multiple tools
process CONTIGS_FROM_IDS {
	tag {"${meta.id}"}
    container 'nakor/virus_extraction'
	publishDir params.outdir, mode: "copy"

	input:
    tuple val(meta), path(fasta), path(ids)

	output:
    tuple val(meta), path("viral-contigs_${meta.id}.fasta")
    
	script:
    def ids_all_tools = ids.join(' ')
    """
    python contigs_from_ids.py \ 
        --fasta $fasta \ 
        --ids $ids_all_tools \
        --output viral-contigs_${meta.id}.fasta
    """
}

// concatenate and clusters the assemblies from multiple samples
// need to add a prefix (e.g. sample name) in case the contig ids
// are the same in the different assemblies
process MERGE_CONTIGS {
    publishDir "${params.outdir}/cluster_genomes", mode: 'copy'
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.fna'
    container 'nakor/votu_table'

    input:
    path f
    val samples

    output:
    path '*', emit: all
    path '*.fna', emit: fasta
    path '*.clstr', emit: clusters    

    when:
    f.size() > 1

    script:
    """
    #!/usr/bin/env bash

    bash add_prefix.sh ${samples.join(' ')}
    cat *_prefixed.fasta > all_contigs.fasta && rm -f *prefixed.fasta
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
    samtools index ${meta}.bam ${meta}.bai
    """
}


process MAKE_COVERAGE_TABLE {
    publishDir "${params.outdir}", mode: 'copy'
    container 'nakor/votu_table'

    input:
    path bam
    path bai

    output:
    path "coverage_table.tsv"

    script:
    """
    #!/usr/bin/env bash

    samtools depth ${bam.sort().join(' ')} > coverage_table.tsv
    """
}

process MAKE_ABUNDANCE_TABLE {
    publishDir "${params.outdir}", mode: 'copy'
    container 'nakor/votu_table'

    input:
    val samples
    path table
    path fasta

    output:
    path "abundance_table.csv"

    script:
    """
    python compute_table.py \
        --fasta ${fasta} \
        --coverage ${table} \
        --samples ${samples.join(' ')}
    """
}


process SUBSET_ABUNDANCE_TABLE {
    publishDir "${params.outdir}", mode: 'copy'
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
    for tool in ${tools_list.join(' ')}; do
        python subset_abundance_table.py \
            --name \$tool \
            --ids *\${tool}*.txt \
            --clusters ${clusters} \
            --abundance ${abundance_table}
    done
    """
}


workflow v_otu_table{
    take:
    sample_names
    contigs
    vir_ids
    reads

    main:

    // subset viral contigs
    // concatenates contigs from all samples
    viral_assembly = CONTIGS_FROM_IDS(
        contigs.combine(vir_ids, by: 0)
    )

    // merge assemblies from multiple samples
    merged_assembly = MERGE_CONTIGS(
        viral_assembly.collect{it[1]}, //
        sample_names
    )

    // align reads against clustered assembly
    ref = BWA_INDEX(merged_assembly.fasta)
    alignments = BWA_MEM(reads, ref)

    // compute coverage
    coverage_table = MAKE_COVERAGE_TABLE(
        alignments.bam.collect(),
        alignments.bai.collect()
    )

    // into abundance table
    abundance_table = MAKE_ABUNDANCE_TABLE(
        sample_names,
        coverage_table,
        merged_assembly.fasta
    )

    // abundance table for each virus identification software
    SUBSET_ABUNDANCE_TABLE(
        vir_ids.collect{it[0]}.unique(),
        vir_ids.collect{it[1]},
        abundance_table,
        merged_assembly.clusters
    )
}

workflow {
    // contig id must be the same as virus ids files
    contigs = Channel.fromPath(params.fasta).map{[[id: "${it.getSimpleName()}"], it]}
    // vir_ids filename are {tool}_contigs_{assembly_prefix}.txt
    vir_ids = Channel.fromPath(params.ids)
        .map{[[id: "${it.getSimpleName().tokenize('_')[2..-1].join('_')}"], it]}
        .groupTuple() // emits tuple(sample, files)
    reads = Channel.fromFilePairs(params.fastq)
    
    sample_names = reads.collect{it[0]}.sort()
    v_otu_table(contigs, vir_ids, reads, sample_names)
}
