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
printf 'Accession\tTolid\tSRA\t# reads\tTotal read length\tAverage read length\tRead N50\tSmallest read length\tLargest read length\tCoverage\tGC content \tBase composition (A:C:T:G)\tAverage read quality\n' >> rdeval.tsv
while IFS="," read -r accession tolid SRA
do
	RANDOM=$SEED
	VAL=$RANDOM
	SEED=$RANDOM
	if (( $(echo "scale=4; ${VAL}/32767 > 0.05" |bc -l) )); then
		printf "skipping: $accession\t$tolid\t$SRA\n"
        continue
    fi
    esearch -db sra -query $SRA | esummary | xtract -pattern DocumentSummary -element Sample@acc Run@acc Experiment@acc Platform instrument_model LIBRARY_STRATEGY Summary -element Statistics@total_bases > accessions.ls
    cat accessions.ls > all_accessions.ls
    
    while read -r SAMPLE SRR SRX INSTR LIB_STR
    do
        printf "downloading: $SAMPLE\t$SRR\t$SRX\n"
        fasterq-dump $SRR
        
        printf "$SAMPLE\t$SRR\t$SRX\t" >> rdeval.tsv
	    rdeval -r ${SRR}*.fastq | cut -d':' -f2 | sed 's/^.//' >> rdeval.tsv
	
	    printf "$SAMPLE\t$SRR\t$SRX\t" >> rdevalCumInv.tsv
	    rdeval ${SRR}*.fastq -s c >> rdevalCumInv.tsv
        
    done<accessions.ls

	rm -f ${SRR}*.fastq
done< <(grep SRS raw_data_metadata.ls)
```