#!/usr/bin/env python3
from Bio import SeqIO
import pandas as pd
import argparse

def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str)
    parser.add_argument('--coverage', type=str)
    parser.add_argument('--samples', type=str, nargs='+')        
    args = parser.parse_args()

    return args

def compute_table(fasta, coverage, sample_names):
    ctg_len = pd.Series(
        {ctg.id: len(ctg.seq) for ctg in SeqIO.parse(fasta, 'fasta')}
    )

    table = pd.read_csv(
        coverage, sep='\t', names=["contig", "pos"] + sample_names
    ).drop(columns="pos")

    table = table.groupby('contig').sum()
    table = (table.T / ctg_len).T.fillna(0)

    table.to_csv('abundance_table.csv')

if __name__ == '__main__':
    args = parse_args()
    compute_table(args.fasta, args.coverage, args.samples)
