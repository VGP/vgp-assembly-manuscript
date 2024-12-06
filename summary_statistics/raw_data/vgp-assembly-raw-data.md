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
We can also run [rdeval](https://github.com/vgl-hub/rdeval) on each data set on the fly to get the summary statistics using [rdeval.sh](rdeval.sh).