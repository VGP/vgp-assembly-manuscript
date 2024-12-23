#!/bin/bash
set -e
. env_parallel.bash

function parallel_download() {
  accession="$1"
  printf "prefetch accession $accession\n"
  counter=1
  while [[ $counter -le 10 ]] ; do
    printf "attempt $counter\n"
    prefetch $accession --max-size u && break
    ((counter++))
  done

  printf "sra to fastq for accession $accession\n"
  counter=1
  while [[ $counter -le 10 ]] ; do
    printf "attempt $counter\n"
    fasterq-dump $accession && break
    ((counter++))
  done
}
export -f parallel_download
rm -f all_accessions.ls
SEED=42
if [ ! -s rdeval.tsv ]; then # add header
  printf 'SRS\t# reads\tTotal read length\tAverage read length\tRead N50\tSmallest read length\tLargest read length\tCoverage\tGC content\tBase composition (A:C:T:G)\tAverage read quality\n' >> rdeval.tsv
fi
while IFS="," read -r -u 3 accession tolid SRA
do
	RANDOM=$SEED
	VAL=$RANDOM
	SEED=$RANDOM
	printf "Processing: %s\t%s\t%s\n" "$accession" "$tolid" "$SRA"
	if (( $(echo "scale=4; ${VAL}/32767 > 0.05" |bc -l) )); then
		printf "Skipping for subsampling.\n"
    continue
  fi
  printf "Searching: %s\n" "$SRA"
  esearch -db sra -query $SRA | esummary | xtract -pattern DocumentSummary -element Sample@acc Run@acc Experiment@acc Platform instrument_model LIBRARY_STRATEGY Summary -element Statistics@total_bases | grep SMRT | grep 'WGS\|WGA' | awk '{if ($6!=0) print}' > accessions.ls
  cat accessions.ls >> all_accessions.ls
  printf "Found records:\n"
  cat accessions.ls
  if grep -q "$SRA" rdeval.tsv; then
		printf "Already done. Skipping.\n"
    continue
  fi

  cat accessions.ls | env_parallel -j 32 --colsep '\t' parallel_download {2}

  printf "Computing summary statistics...\n"
  printf "%s\t" "$SRA" >> rdeval.tsv
  rdeval -r *.fastq | awk -F': ' '{print $2}' | sed 1d | sed -z 's/\n/\t/g; s/.$//' >> rdeval.tsv
  printf "\n" >> rdeval.tsv

  printf "Computing Cumulative inverse distribution...\n"
  printf "%s\t" "$SRA" >> rdevalCumInv.tsv
  rdeval *.fastq -s c | sed -z 's/\n/;/g' >> rdevalCumInv.tsv
  printf "\n" >> rdevalCumInv.tsv

	rm -f *.fastq
done 3< <(grep 'ERS\|SRS' raw_data_metadata.ls | grep -v alt)