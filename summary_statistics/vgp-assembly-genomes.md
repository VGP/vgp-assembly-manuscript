# VGP assembly exploratory analysis
NCBI's dataset utility can be downloaded from [here](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).

Download VGP Bioproject metadata:

```
datasets summary genome accession PRJNA489243 > vgp-metadata.json
```

## list VGP repo

We can use [jq](https://jqlang.github.io/jq/) to parse the repo (here is jq's [manual](https://jqlang.github.io/jq/manual/)):
```
cat vgp-metadata.json | jq -r '.reports[] | .accession + "," + .assembly_info.assembly_name + "," + (.assembly_info.biosample.attributes[] | select(.name=="scientific_name").value) // .accession + "," + .assembly_info.assembly_name + ","' > accession_metadata.ls
cat accession_metadata.ls
```

## Download assemblies

We can download a random subset of genomes combining jq's and NCBI's datasets functionalities.
We can also run [gfastats](https://github.com/vgl-hub/gfastats) on each genome on the fly to get the summary statistics.

```
#!/bin/bash
set -e

SEED=42
rm -f gfastats.tsv gfastatsNxContig.tsv gfastatsNxScaffold.tsv
printf 'Accession\tTolid\tScientific name\t# scaffolds\tTotal scaffold length\tAverage scaffold length\tScaffold N50\tScaffold auN\tScaffold L50\tLargest scaffold\tSmallest scaffold\t# contigs\tTotal contig length\tAverage contig length\tContig N50\tContig auN\tContig L50\tLargest contig\tSmallest contig\t# gaps in scaffolds\tTotal gap length in scaffolds\tAverage gap length in scaffolds\tGap N50 in scaffolds\tGap auN in scaffolds\tGap L50 in scaffolds\tLargest gap in scaffolds\tSmallest gap in scaffolds\tBase composition (A\tGC content\t# soft-masked bases\t# segments\tTotal segment length\tAverage segment length\t# gaps\t# paths\n' >> gfastats.tsv
while IFS="," read -r accession tolid latin_name
do
	RANDOM=$SEED
	VAL=$RANDOM
	SEED=$RANDOM
	if (( $(echo "scale=4; ${VAL}/32767 > 0.25" |bc -l) )); then
		printf "skipping: $accession\t$tolid\t$latin_name\n"
        continue
    fi
	datasets download genome accession $accession --filename $accession.zip
	unzip -o $accession.zip
	printf "$accession\t$tolid\t$latin_name\t" >> gfastats.tsv
	
	genome=ncbi_dataset/data/$accession/*.fna
	if ! cmp -s <(md5sum $genome) <(grep fna md5sum.txt); then
		printf "Check file integrity: $accession"
		exit 1
	fi
	
	gfastats -t $genome | cut -f2 | sed -z 's/\n/\t/g; s/.$//' >> gfastats.tsv
	
	printf "$accession\t$tolid\t$latin_name\t" >> gfastatsNxContig.tsv
	gfastats $genome -s c | sort -nrk2 | awk 'BEGIN{pos=0}{total+=$2; size[pos] = $2; cum_size[pos++] = total}END{for (p = 0; p < pos; p++) {printf size[p]","cum_size[p]/total"\t"}; printf "\n"}' >> gfastatsNxContig.tsv
	
	printf "$accession\t$tolid\t$latin_name\t" >> gfastatsNxScaffold.tsv
	gfastats $genome -s s | sort -nrk2 | awk 'BEGIN{pos=0}{total+=$2; size[pos] = $2; cum_size[pos++] = total}END{for (p = 0; p < pos; p++) {printf size[p]","cum_size[p]/total"\t"}; printf "\n"}' >> gfastatsNxScaffold.tsv
	
	rm -r $accession.zip ncbi_dataset README.md md5sum.txt
done<accession_metadata.ls
```
