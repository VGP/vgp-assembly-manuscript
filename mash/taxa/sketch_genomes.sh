#!/bin/bash
set -e

SUBSET=${1:-false}
THREAD=${2:-32}

SEED=42
mkdir -p sketches
rm -f genome_list.tsv mash_distances.tsv
while IFS="," read -r accession tolid latin_name
do
	RANDOM=$SEED
	VAL=$RANDOM
	SEED=$RANDOM
	if [ $SUBSET ] && (( $(echo "scale=4; ${VAL}/32767 > 0.25" |bc -l) )); then
		printf "Subsampling mode. Skipping: %s\t%s\t%s\n" "$accession" "$tolid" "$latin_name"
    continue
  fi
  if [ -f sketches/$accession*.msh ]; then
		printf "Sketch already available. Skipping: %s\t%s\t%s\n" "$accession" "$tolid" "$latin_name"
    continue
  fi
	datasets download genome accession $accession --filename $accession.zip
	unzip -o $accession.zip -d $accession
	printf "%s\t%s\t%s\n" "$accession" "$tolid" "$latin_name" >> genome_list.tsv

	genome="$accession/ncbi_dataset/data/$accession/*.fna"
	if ! cmp -s <(md5sum $genome | cut -f1 -d' ') <(grep fna $accession/md5sum.txt | cut -f1 -d' '); then
		printf "Check file integrity: %s" "$accession"
		exit 1
	fi

	mash sketch -p $THREAD -s 10000000 $genome
	mv $genome.msh sketches
	rm -r $accession.zip $accession
done<accession_metadata.ls