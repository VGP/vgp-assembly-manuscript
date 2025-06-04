# VGP assembly exploratory analysis

NCBI's dataset utility can be downloaded from [here](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).

First we retrieve the information for VGP Phase I genomes (primary only):
```
wget https://raw.githubusercontent.com/VGP/vgp-phase1/refs/heads/main/VGPPhase1-freeze-1.0.tsv
awk -F'\t' 'NR>1 {if($16 != "")printf $16","$14","$10","$2"\n"}' VGPPhase1-freeze-1.0.tsv > accession_metadata.ls
```

Alternatively, we can use [jq](https://jqlang.github.io/jq/) to parse the repo (here is jq's [manual](https://jqlang.github.io/jq/manual/)) and download VGP Bioproject metadata:

```
datasets summary genome accession PRJNA489243 > vgp-metadata.json
cat vgp-metadata.json | jq -r '.reports[] | .accession + "," + .assembly_info.assembly_name + "," + (.assembly_info.biosample.attributes[] | select(.name=="scientific_name").value) // .accession + "," + .assembly_info.assembly_name + ","' > accession_metadata.ls
cat accession_metadata.ls
```

## Download assemblies

We can download a random subset of genomes combining jq's and NCBI's datasets functionalities (default: download all, SUBSAMPLING_FRACTION).
We can also run [gfastats](https://github.com/vgl-hub/gfastats) on each genome on the fly to get the summary statistics. [gfastats.sh](gfastats.sh) does that, including also Nx plots.

## Plotting
There are a number of scripts to plot:
- [Nx_plots.py](Nx_plots.py) for Nx statistics
