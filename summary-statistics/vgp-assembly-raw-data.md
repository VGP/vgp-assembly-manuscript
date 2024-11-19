# VGP assembly exploratory analysis
NCBI's dataset utility can be downloaded from [here](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).

Download VGP Bioproject metadata:

```
datasets summary genome accession PRJNA489243 > vgp-metadata.json
```

## list VGP repo

We can use [jq](https://jqlang.github.io/jq/) to parse the repo (here is jq's [manual](https://jqlang.github.io/jq/manual/)):
```
cat vgp-metadata.json | jq -r '.reports[] | .accession + "," + .assembly_info.assembly_name + "," + (.assembly_info.biosample.sample_ids[] | select(.db=="SRA").value)' > raw_data_metadata.ls
cat raw_data_metadata.ls
```

## Download raw data

We can download a random subset of raw data sets combining jq's and NCBI's Entrez Direct ([EDirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/)) functionalities.
We can also run [rdeval](https://github.com/vgl-hub/rdeval) on each data set on the fly to get the summary statistics.

```
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
	printf "Processing: $accession\t$tolid\t$SRA\n"
	if (( $(echo "scale=4; ${VAL}/32767 > 0.05" |bc -l) )); then
		printf "Skipping for subsampling.\n"
        continue
    fi
    esearch -db sra -query $SRA | esummary | xtract -pattern DocumentSummary -element Sample@acc Run@acc Experiment@acc Platform instrument_model LIBRARY_STRATEGY Summary -element Statistics@total_bases > accessions.ls
    cat accessions.ls >> all_accessions.ls
    
    while read -r -u 4 SAMPLE SRR SRX INSTR LIB_STR
    	do
        	printf "downloading: $SAMPLE\t$SRR\t$SRX\n"
        	fasterq-dump $SRR
        done 4< <(grep SMRT accessions.ls)
        
    printf "$SAMPLE\t" >> rdeval.tsv
    rdeval -r *.fastq | awk -F': ' '{print $2}' | sed 1d | sed -z 's/\n/\t/g; s/.$//' >> rdeval.tsv

    printf "$SAMPLE\t" >> rdevalCumInv.tsv
    rdeval *.fastq -s c | sed -z 's/\n/;/g' >> rdevalCumInv.tsv

	rm -f *.fastq
done 3< <(grep 'ERS\|SRS' raw_data_metadata.ls)
```