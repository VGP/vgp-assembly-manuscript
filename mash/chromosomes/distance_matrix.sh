#!/bin/bash

THREAD="${1:-32}"

rm -f distance_matrix.tsv

readarray -t in < $1 # input list

for (( i=0; i<${#in[@]}; i++ )); do # triangular matrix
    for (( j=0; j<${#in[@]}; j++ )); do
        if [ $j -gt $i ]; then # sketch and get the dist
        	accession1=${in[$i]%%$'\t'*}
        	accession2=${in[$j]%%$'\t'*}
          awk -v acc1="$accession1" -v acc2="$accession2" '{printf acc1"\t"acc2"\t"$0"\n"}' <(mash dist -p "$THREAD" sketches/$accession1* sketches/$accession2*) >> distance_matrix.tsv
        fi
    done
done