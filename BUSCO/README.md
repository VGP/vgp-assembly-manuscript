# BUSCO analysis

Contained in this folder are statistics from BUSCO runs as well as contiguity statistics for the VGP phase I genomes as well as genomes published on ncbi under the DNAZoo and Zoonomia (200 mammals) umbrella projects.

BUSCO was run in genome mode using the vertebrata_odb10 database and miniprot mapper like so:

```apptainer exec -B /scratch/brown/ /scratch/brown/progs/busco_v5.7.1_cv1.sif busco -m genome -l vertebrata --download_path /path/th/db/ --offline -o ncbi_dataset/data/${f}/${f}.miniprot.busco -c 12 -i ncbi_dataset/data/${f}/*fna```

read `singularity` in place of `apptainer` where appropriate.

The script `get_stats.sh` is also used to extract BUSCO, n50, size, number statistics as well as look up taxonomic information about each species. The conda environment `jq.yaml` includes the software used to query the xml files coming from the ena browser.
