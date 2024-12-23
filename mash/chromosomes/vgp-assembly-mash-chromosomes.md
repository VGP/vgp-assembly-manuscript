# VGP assembly mash graph
## Download VGP metadata from NCBI
NCBI's dataset utility can be downloaded from [here](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).
Download VGP Bioproject metadata:
```
datasets summary genome accession PRJNA489243 > vgp-metadata.json
```
## list VGP repo
We can use [jq](https://jqlang.github.io/jq/) to parse the repo (here is jq's [manual](https://jqlang.github.io/jq/manual/)), selecting only primary/complete genomes:
```
cat vgp-metadata.json | 
jq -r '.reports[] | 
.accession + "," + .assembly_info.assembly_name + "," + '.assembly_info.assembly_type' + "," + .assembly_info.diploid_role + "," + .assembly_info.assembly_level + "," + (.assembly_info.biosample.attributes[] | select(.name=="scientific_name").value) // 
.accession + "," + .assembly_info.assembly_name + "," + '.assembly_info.assembly_type' + "," + .assembly_info.diploid_role + "," + .assembly_info.assembly_level' | \
grep -v alt | grep Chromosome > accession_metadata.ls
cat accession_metadata.ls
```
Note: we used jq's alternative operator `\\` in case species not set
## Download assemblies
We can download a random subset of genomes combining jq's and NCBI's datasets functionalities, then compute [mash](https://github.com/marbl/Mash) sketches and all-vs-all distances.
First compute individual mash sketches:
```
#!/bin/bash
set -e

SEED=42
mkdir -p sketches
rm -f genome_list.tsv mash_distances.tsv
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
```
Next compute triangular mash distance matrix with [distance_matrix.sh](../distance_matrix.sh):
```
bash distance_matrix.sh genome_list.tsv
```

We can now run `distance_hist.py` to get the histogram of the distances.
We can add metadata to the accessions using this command:
```
awk 'FNR==NR{a[$1]=$2; next} {FS=" "} {$0=$0; print $1","a[$1]","$2","a[$2]","$3","$4","$5}' FS="," accession_metadata.ls distance_matrix.tsv > distance_matrix_with_accessions.csv
```
Now we can find individual chromosome matches with [best_hits.py](best_hits.py)