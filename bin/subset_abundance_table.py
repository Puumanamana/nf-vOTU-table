#!/usr/bin/env python

from pathlib import Path
import pandas as pd
import argparse

def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--name', type=str, default='tool')
    parser.add_argument('--abundance', type=str)
    parser.add_argument('--clusters', type=str, default='')
    parser.add_argument('--ids', type=str, nargs='+')        
    args = parser.parse_args()

    return args

def parse_clstr_file(clstr_file):
    
    clusters = {}
    repr_seq = None

    with open(clstr_file, 'r') as handle:
        for line in handle:
            entries = line.strip().split('\t')
            if line.startswith('>'):
                repr_seq = entries[1]
                clusters[repr_seq] = repr_seq
            else:
                clusters[entries[0]] = repr_seq

    return pd.Series(clusters)

if __name__ == '__main__':
    args = parse_args()

    # Write the abundace tables
    abd = pd.read_csv(args.abundance, index_col=0)

    ids = {'{}!{}'.format(
        Path(f).stem.split('-')[-1],
        x.strip().replace('.', '_').split('\t')[0]
    ) for f in args.ids for x in open(f)}

    if Path(args.clusters).is_file():
        clstrs = parse_clstr_file(args.clusters)
        clstrs.index = clstrs.index.str.replace('.', '_') # for Virsorter
        votu_subset = clstrs[list(ids)].unique()
    else:
        votu_subset = [x.split()[0] for x in ids]

    abd.loc[votu_subset].to_csv(f'abundance_table-{args.name}.csv') 
