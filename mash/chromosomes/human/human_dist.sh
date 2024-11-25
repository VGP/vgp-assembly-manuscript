#!/bin/bash

THREAD="${2:-32}"

rm -f human_distances.tsv

readarray -t in < $1 # input list

for (( i=0; i<${#in[@]}; i++ )); do
	accession=${in[$i]%%$'\t'*}
	awk -v acc="$accession" '{printf acc"\t"$0"\n"}' <(mash dist -p "$THREAD" $2 $3/$accession*) >> human_distances.tsv
done
