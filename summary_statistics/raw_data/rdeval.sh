#!/bin/bash
set -e

SEED=42
rm -f rdeval.tsv rdevalCumInv.tsv all_accessions.ls
printf 'SRS\tSRR\tSRX\t# reads\tTotal read length\tAverage read length\tRead N50\tSmallest read length\tLargest read length\tCoverage\tGC content\tBase composition (A:C:T:G)\tAverage read quality\n' >> rdeval.tsv
while IFS="," read -r -u 3 accession tolid SRA
do
	RANDOM=$SEED
	VAL=$RANDOM
	SEED=$RANDOM
	printf "Processing: %s\t%s\t%s\n", "$accession", "$tolid", "$SRA"
	if (( $(echo "scale=4; ${VAL}/32767 > 0.05" |bc -l) )); then
		printf "Skipping for subsampling.\n"
    continue
  fi
  printf "Searching: %s\n", "$SRA"
  esearch -db sra -query "$SRA" | esummary | xtract -pattern DocumentSummary -element Sample@acc Run@acc Experiment@acc Platform instrument_model LIBRARY_STRATEGY Summary -element Statistics@total_bases > accessions.ls
  cat accessions.ls >> all_accessions.ls
  printf "Found records:\n"
  cat accessions.ls

  while read -r -u 4 SAMPLE SRR SRX INSTR LIB_STR
    do
        if [ ! -f $SRR.fastq ]; then
          printf "downloading: $SAMPLE\t$SRR\t$SRX\n", "$SAMPLE", "$SRR", "$SRX"
          fasterq-dump $SRR
        fi
    done 4< <(grep SMRT accessions.ls | awk '{if ($6!=0) print}')

  printf "$SAMPLE\t" >> rdeval.tsv
  rdeval -r *.fastq | awk -F': ' '{print $2}' | sed 1d | sed -z 's/\n/\t/g; s/.$//' >> rdeval.tsv

  printf "$SAMPLE\t" >> rdevalCumInv.tsv
  rdeval *.fastq -s c | sed -z 's/\n/;/g' >> rdevalCumInv.tsv

	rm -f *.fastq
done 3< <(grep 'ERS\|SRS' raw_data_metadata.ls)