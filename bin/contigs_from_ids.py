from Bio import SeqIO
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str)
    parser.add_argument('--ids', type=str, nargs='+')
    parser.add_argument('--output', type=str)    
    args = parser.parse_args()

    return args

def extract_contigs(fasta, output, *txt):
    viral_ids = {ctg.strip().split(',')[0].split()[0] for f in txt for ctg in open(f)}

    viral_sequences = []
    for ctg in SeqIO.parse(fasta, 'fasta'):
        if ctg.id in viral_ids:
            viral_sequences.append(ctg)

    SeqIO.write(viral_sequences, output, 'fasta-2line')

if __name__ == '__main__':
    args = parse_args()
    extract_contigs(args.fasta, args.output, args.ids)
