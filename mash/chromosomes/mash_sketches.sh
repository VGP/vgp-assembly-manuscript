#!/bin/bash
# requires datasets, gfastata, unzip, mash, md5sum
set -e

SUBSAMPLE=$1 # subsbampling fraction
printf "Using subsbampling fraction: ${SUBSAMPLE}\n"

SEED=42
mkdir -p sketches
rm -f genome_list.tsv mash_distances.tsv
while IFS="," read -r accession tolid latin_name
do
	RANDOM=$SEED
	VAL=$RANDOM
	SEED=$RANDOM
	if (( $(echo "scale=4; ${VAL}/32767 > ${SUBSAMPLE}" |bc -l) )); then
		printf "Skipping: $accession\t$tolid\t$latin_name\n"
        continue
    fi
	datasets download genome accession $accession --filename $accession.zip
	unzip -o $accession.zip -d $accession
	printf "$accession\t$tolid\t$latin_name\n" >> genome_list.tsv

	genome=$accession/ncbi_dataset/data/$accession/*.fna
	if ! cmp -s <(md5sum $genome | cut -f1 -d' ') <(grep fna $accession/md5sum.txt | cut -f1 -d' '); then
		printf "Check file integrity: $accession"
		exit 1
	fi

	gfastats -s s $genome | cut -f1 > chr.ls # assembled molecules only
	gfastats -i <(grep $(cut -c1-2 chr.ls | uniq | head -1) chr.ls) $genome -o $accession.chr.fa
	mash sketch -p 32 -i -s 10000000 $accession.chr.fa # parellized 32 ways, per chromosome, 10M hashes
	mv $accession.chr.fa.msh sketches
	rm -r $accession.chr.fa $accession.zip $accession
done<accession_metadata.ls