#!/usr/bin/env python

from Bio import SeqIO
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str)
    parser.add_argument('--ids', type=str, nargs='+')
    parser.add_argument('--output', type=str)    
    parser.add_argument('--prefix', type=str, default='')    
    args = parser.parse_args()

    return args

def extract_contigs(fasta, output, ids=None, prefix=''):
    # We need to replace '.' in names because virsorter modifies ctg ids
    viral_ids = {ctg.strip().split(',')[0].split()[0].replace('.', '_')
                 for f in ids for ctg in open(f)}

    viral_sequences = []
    for ctg in SeqIO.parse(fasta, 'fasta'):
        if ctg.id.replace('.', '_') in viral_ids:
            ctg.id = f'{prefix}{ctg.id}'
            ctg.description = ''
            viral_sequences.append(ctg)

    SeqIO.write(viral_sequences, output, 'fasta-2line')

if __name__ == '__main__':
    args = parse_args()
    extract_contigs(args.fasta, args.output, ids=args.ids, prefix=args.prefix)
