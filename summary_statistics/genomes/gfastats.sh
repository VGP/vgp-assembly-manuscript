#!/bin/bash
set -e
. env_parallel.bash

SUBSAMPLE=$1 # subsampling fraction

function parallel_gfastats() {
  accession="$1"
  tolid="$2"
  latin_name="$3"

  # Create a subfolder for each accession
  mkdir -p accessions/$accession
  cd accessions/$accession

  # Download and unzip the genome into the accession-specific folder
  datasets download genome accession $accession --filename $accession.zip
  unzip -o $accession.zip

  # Generate individual files for each accession
  printf "$accession\t$tolid\t$latin_name\t" > gfastats.tsv
  gfastats -t ncbi_dataset/data/$accession/*.fna | cut -f2 | sed -z 's/\n/\t/g; s/.$//' >> gfastats_$accession.tsv
  echo "" >> gfastats_$accession.tsv
  printf "$accession\t$tolid\t$latin_name\t" > gfastatsNxContig.tsv
  gfastats ncbi_dataset/data/$accession/*.fna -s c | sort -nrk2 | awk 'BEGIN{pos=0}{total+=$2; size[pos] = $2; cum_size[pos++] = total}END{for (p = 0; p < pos; p++) {printf size[p]","cum_size[p]/total"\t"}; printf "\n"}' >> gfastatsNxContig_$accession.tsv
  printf "$accession\t$tolid\t$latin_name\t" > gfastatsNxScaffold_$accession.tsv
  gfastats ncbi_dataset/data/$accession/*.fna -s s | sort -nrk2 | awk 'BEGIN{pos=0}{total+=$2; size[pos] = $2; cum_size[pos++] = total}END{for (p = 0; p < pos; p++) {printf size[p]","cum_size[p]/total"\t"}; printf "\n"}' >> gfastatsNxScaffold_$accession.tsv

  rm -r $accession.zip ncbi_dataset README.md md5sum.txt
  cd ../..
    printf "Processed accession: $accession\n"
}

export -f parallel_gfastats

# Initialize empty files for final output
rm -f gfastats.tsv gfastatsNxContig.tsv gfastatsNxScaffold.tsv
printf 'Accession\tTolid\tScientific name\t# scaffolds\tTotal scaffold length\tAverage scaffold length\tScaffold N50\tScaffold auN\tScaffold L50\tLargest scaffold\tSmallest scaffold\t# contigs\tTotal contig length\tAverage contig length\tContig N50\tContig auN\tContig L50\tLargest contig\tSmallest contig\t# gaps in scaffolds\tTotal gap length in scaffolds\tAverage gap length in scaffolds\tGap N50 in scaffolds\tGap auN in scaffolds\tGap L50 in scaffolds\tLargest gap in scaffolds\tSmallest gap in scaffolds\tBase composition (A\tGC content\t# soft-masked bases\t# segments\tTotal segment length\tAverage segment length\t# gaps\t# paths\n' >> gfastats.tsv

# Process all accessions and apply random subsampling
cat accession_metadata.ls | env_parallel -j 8 --colsep ',' '
  accession={1}; tolid={2}; latin_name="{3}";

  # Random subsampling logic
  SUBSAMPLING_THRESHOLD=$(echo "$SUBSAMPLE * 32767" | bc)
  RAND_NUM=$((RANDOM % 32767))

  if (( RAND_NUM > SUBSAMPLING_THRESHOLD )); then
    printf "Skipping $accession for subsampling.\n"
    continue
  fi

  # Run the function for each accession
  parallel_gfastats $accession $tolid "$latin_name"
'

# After all parallel tasks finish, combine the individual files into the final output files
cat accessions/*/gfastats_*.tsv >> gfastats.tsv
cat accessions/*/gfastatsNxContig_*.tsv >> gfastatsNxContig.tsv
cat accessions/*/gfastatsNxScaffold_*.tsv >> gfastatsNxScaffold.tsv