#!/bin/bash

rm -f mash_chromosomes.out

readarray -t in < $1 # input list

for (( i=0; i<${#in[@]}; i++ )); do # triangular matrix
    for (( j=0; j<${#in[@]}; j++ )); do
        if [ $j -gt $i ]; then # sketch and get the dist
        	accession1=${in[$i]%%$'\t'*}
        	accession2=${in[$j]%%$'\t'*}
          awk -v acc1="$accession1" -v acc2="$accession2" '{printf acc1"\t"acc2"\t"$0"\n"}' <(mash dist sketches/$accession1* sketches/$accession2*) >> mash_chromosomes.out
        fi
    done
done