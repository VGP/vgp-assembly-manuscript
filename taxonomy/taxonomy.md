# VGP assembly mash graph
## Download VGP metadata from NCBI
NCBI's dataset utility can be downloaded from [here](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).
Download VGP Bioproject metadata:
```
datasets summary genome accession PRJNA489243 > vgp-metadata.json
```
## list VGP repo
We can use [jq](https://jqlang.github.io/jq/) to parse the repo (here is jq's [manual](https://jqlang.github.io/jq/manual/)) and get all accessions with species associated with the Biosample:
```
cat vgp-metadata.json | jq -r '.reports[] | .accession + "," + .assembly_info.assembly_name + "," + (.assembly_info.biosample.attributes[] | select(.name=="scientific_name").value) // .accession + "," + .assembly_info.assembly_name + ","' > accession_metadata.ls
cat accession_metadata.ls
```
Get all species with accessions
```
cat vgp-metadata.json | jq -r '.reports[] | .accession + "," + .assembly_info.assembly_name + "," + (.assembly_info.biosample.attributes[] | select(.name=="scientific_name").value)' | cut -d',' -f3 | sort | uniq > species.ls
cat species.ls
```
We can use R's [taxonomizr](https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html) package to explore the taxonomy further:
```
Rscript taxonomizr.R
```