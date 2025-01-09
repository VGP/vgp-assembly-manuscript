#!/bin/bash
set -e
. env_parallel.bash

SUBSAMPLE=$1 # subsbampling fraction

function parallel_rdeval() {
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

  if [ -s $accession.fastq ]; then # check if file exists
    printf "Generating .rd file and summary statistics...\n"
    printf "%s\t" "$accession" >> rdeval_$accession.tsv
    rdeval $accession.fastq -o $accession.rd | awk -F': ' '{print $2}' | sed 1d | sed -z 's/\n/\t/g; s/.$//' >> rdeval_$accession.tsv
    printf "\n" >> rdeval_$accession.tsv
    printf "Computing Cumulative inverse distribution...\n"
    printf "%s\t" "$accession" >> rdevalCumInv_$accession.tsv
    rdeval $accession.rd -s c | sed -z 's/\n/;/g' >> rdevalCumInv_$accession.tsv
    printf "\n" >> rdevalCumInv_$accession.tsv
  else # record accessions with no records
    printf "%s\t" "$accession" >> missing_$accession
  fi
}
export -f parallel_rdeval
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
	if (( $(echo "scale=4; ${VAL}/32767 > ${SUBSAMPLE}" |bc -l) )); then
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

  mkdir -p $SRA
  cd $SRA
  cat ../accessions.ls | env_parallel -j 8 --colsep '\t' parallel_rdeval {2}
  cd ..

  printf "Combining summary statistics...\n"
  printf "%s\t" "$SRA" >> rdeval.tsv
  rdeval $SRA/*.rd | awk -F': ' '{print $2}' | sed 1d | sed -z 's/\n/\t/g; s/.$//' >> rdeval.tsv
  printf "\n" >> rdeval.tsv

  printf "Combining Cumulative inverse distribution...\n"
  printf "%s\t" "$SRA" >> rdevalCumInv.tsv
  rdeval $SRA/*.rd -s c | sed -z 's/\n/;/g' >> rdevalCumInv.tsv
  printf "\n" >> rdevalCumInv.tsv

	rm -f $SRA/*.fastq
done 3< <(grep 'ERS\|SRS' raw_data_metadata.ls | grep -v alt)