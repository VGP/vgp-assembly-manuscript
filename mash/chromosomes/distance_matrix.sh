#!/bin/bash

readarray -t in < $1 # input list

for (( i=0; i<${#in[@]}; i++ )); do # triangular matrix
    for (( j=0; j<${#in[@]}; j++ )); do
        if [ $j -gt $i ]; then # sketch and get the dist
        	accession1=${in[$i]%%$'\t'*}
        	accession2=${in[$j]%%$'\t'*}
          awk -v acc1="$accession1" -v acc2="$accession2" '{print $acc1 $acc2 $0}' <(mash dist sketches/$accession1* sketches/$accession2*)
        fi
    done
done