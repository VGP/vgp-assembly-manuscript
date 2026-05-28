# VGP assembly mash graph
First we retrieve the information for VGP Phase I genomes (primary only):
```
wget https://raw.githubusercontent.com/VGP/vgp-phase1/refs/heads/main/VGPPhase1-freeze-1.0.tsv
awk -F'\t' 'NR>1 {if($16 != "")printf $16","$14","$10","$2"\n"}' VGPPhase1-freeze-1.0.tsv > accession_metadata.csv
```
Alternatively we can use NCBI's VGP project (see below).

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
We can download a random subset of genomes combining jq's and NCBI's datasets functionalities, then compute [mash](https://github.com/marbl/Mash) sketches and all-vs-all distances with [mash_sketches.sh](mash_sketches.sh).
```
bash mash_sketches.sh accession_metadata.csv 1 # subsampling fraction
```

Next compute triangular mash distance matrix with [distance_matrix.sh](../distance_matrix.sh):
```
bash distance_matrix.sh accession_metadata.csv
```

We can now run `distance_hist.py` to get the histogram of the distances.
We can add metadata to the accessions using this command:
```
awk 'FNR==NR{a[$1]=$2; next} {FS=" "} {$0=$0; print $1","a[$1]","$2","a[$2]","$3","$4","$5}' FS="," accession_metadata.ls distance_matrix.tsv > distance_matrix_with_accessions.csv
```
Now we can find individual chromosome matches with [best_hits.py](best_hits.py)