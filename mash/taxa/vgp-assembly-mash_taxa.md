# VGP assembly mash graph
## Download VGP metadata from NCBI
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
We can download a random subset of genomes combining jq's and NCBI's datasets functionalities, then compute [mash](https://github.com/marbl/Mash) sketches and all-vs-all distances.
First compute individual mash sketches with [sketch_genomes.sh](sketch_genomes.sh):
```
bash sketch_genomes.sh false 32 # do not subset, nthreads
```
Next compute triangular mash distance matrix with [distance_matrix.sh](distance_matrix.sh):
```
bash distance_matrix.sh genome_list.tsv
```

We can now run `distance_hist.py` to get the histogram of the distances.
We can add metadata to the accessions using this command:
```
awk 'FNR==NR{a[$1]=$2; next} {FS=" "} {$0=$0; print $1, a[$1], a[$2], $2, $3}' FS="," accession_metadata.ls distance_matrix.txt
```
Potentially filtering distant interactions:
```
printf "Accession 1,Tolid 1,Class 1,Accession 2,Tolid 2,Class 2,D\n" > filtered_distance_matrix.txt
awk 'FNR==NR{a[$1]=$2; next} {FS=" "} {$0=$0; if ($3<0.2) printf $1","a[$1]","substr(a[$1], 1, 1)","$2","a[$2]","substr(a[$2], 1, 1)","$3"\n"}' FS="," accession_metadata.ls distance_matrix.txt >> filtered_distance_matrix.txt
```
This can then be visualized in tools such as Cytoscape.
We can also generate a heatmap/hierarchical clustering using the distance matrix:
```
python distance_heatmap.py # work in progress
```