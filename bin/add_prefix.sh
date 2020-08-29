#!/usr/bin/env bash

for s in $@; do
    ids=$(ls *${s}*.fasta)

    if [ -z $ids ]; then
        echo 'Error: some assembly files do not contain the sample name in their name'
        return 1
    elif [ $ids == *' '* ]; then
         echo 'Error: some assembly files contain the same sample name in their name'
         return 1
    fi

    sed -e 's/^/${s}!/' $ids > ${s}_prefixed.fasta
done
