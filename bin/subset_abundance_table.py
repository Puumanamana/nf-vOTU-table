import pandas as pd
import argparse

def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--name', type=str, default='tool')
    parser.add_argument('--abundance', type=str)
    parser.add_argument('--clusters', type=str)
    parser.add_argument('--ids', type=str, nargs='+')        
    args = parser.parse_args()

    return args

def parse_clstr_file(clstr_file):
    
    clusters = {}
    repr_seq = None

    with open(clstr_file, 'r') as handle:
        for line in handle:
            entries = line.split('\t')
            if line.startswith('>'):
                repr_seq = entries[1]
                clusters[repr_seq] = repr_seq
            else:
                clusters[repr_seq] = entries[0]

    return pd.Series(clusters)

if __name__ == '__main__':
    args = parse_args()
    clstrs = parse_clstr_file(args.clusters)

    # Write the abundace tables
    abd = pd.read_csv(args.abundance, index_col=0)

    ids = {f'{args.name}!{x}' for f in args.ids for x in open(f)}
    abd.loc[clstrs[ids]].to_csv(f'abundance-table_{args.name}.csv') 
