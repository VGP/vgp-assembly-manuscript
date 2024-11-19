#!/bin/bash

readarray -t in < $1 # input list

for (( i=0; i<${#in[@]}; i++ )); do # triangular matrix
    for (( j=0; j<${#in[@]}; j++ )); do
        if [ $j -gt $i ]; then # sketch and get the dist
        	accession1=${in[$i]%%$'\t'*}
        	accession2=${in[$j]%%$'\t'*}
            echo $accession1 $accession2 $(mash dist sketches/$accession1* sketches/$accession2* | cut -f3-5)
        fi
    done
done