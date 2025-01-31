#!/bin/bash

THREAD="${2:-32}"

readarray -t in < $1 # input list

for (( i=0; i<${#in[@]}; i++ )); do # triangular matrix
    for (( j=0; j<${#in[@]}; j++ )); do
        if [ $j -gt $i ]; then # sketch and get the dist
        	accession1=${in[$i]%%$'\t'*}
        	accession2=${in[$j]%%$'\t'*}

        	printf "Processing: ${accession1}\t${accession2}\n"
          if grep -qwE "${accession1}.*${accession2}|${accession2}.*${accession1}" distance_matrix.tsv; then
            printf "Already done. Skipping.\n"
            continue
          fi
          awk -v acc1="$accession1" -v acc2="$accession2" '{printf acc1"\t"acc2"\t"$0"\n"}' <(mash dist -p "$THREAD" sketches/$accession1* sketches/$accession2*) >> distance_matrix.tsv
        fi
    done
done