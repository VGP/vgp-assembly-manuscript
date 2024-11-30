#!/bin/bash
set -e

function parallel_download() {
  accession="$1"
   printf "downloading accession $accession\n"
  prefetch $accession --max-size u
  fasterq-dump $accession
}
export -f parallel_download

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
  esearch -db sra -query $SRA | esummary | xtract -pattern DocumentSummary -element Sample@acc Run@acc Experiment@acc Platform instrument_model LIBRARY_STRATEGY Summary -element Statistics@total_bases > accessions.ls
  cat accessions.ls >> all_accessions.ls
  printf "Found records:\n"
  cat accessions.ls
  if grep -q "$SRA" rdeval.tsv; then
		printf "Already done. Skipping.\n"
    continue
  fi

  grep SMRT accessions.ls | grep WGS accessions.ls | awk '{if ($6!=0) print}' | parallel -j 32 --colsep '\t' parallel_download ::: {2}

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