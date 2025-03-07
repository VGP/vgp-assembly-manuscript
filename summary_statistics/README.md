# Summary statistics
Files in this folder:
- [vgp-assembly-genomes.md](vgp-assembly-genomes.md) showcases how to access genome accession numbers and metadata from NCBI, download the assembly an generate summary statistics
- [vgp-assembly-raw-data.md](vgp-assembly-raw-data.md) showcases how to access raw data associated with individual genomes and generate summary statistics
- [vgp-assembly-genomes.py](vgp-assembly-genomes.py) generates analyses and figures using the summary statistics on genomes

Handy commands to extract all non-VGP genomes with certain metrics:

```
datasets summary genome taxon 7742 > all_vertebrates.07032025.json
cat all_vertebrates.07032025.json | jq -r '.reports[] | select(all(.assembly_info.bioproject_lineage[].bioprojects[].parent_accessions[]?; . != "PRJNA489243") and .assembly_info.assembly_level == "Chromosome" and .assembly_stats.contig_n50 > 1000000 and .assembly_stats.scaffold_n50 > 10000000 and (.assembly_info.assembly_method | index("hifiasm"))) | [.accession, .assembly_info.assembly_name, .assembly_info.assembly_level, (.assembly_stats.contig_n50|tostring), (.assembly_stats.scaffold_n50|tostring), .assembly_info.assembly_method] | join(",")'
```
