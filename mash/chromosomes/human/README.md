# label human ortholog chromosomes
To label human chromosomes use `human_dist.sh` to compute pairwise distances between human chromosome sketches and candidate chromosomes:
```
bash human_dist.sh accessions.ls human.fasta.msh sketches/
```
Note: `sketches/` is the folder where all the sketches that need to be labelled are, prefixed with the accession number.
It generates `human_distances.tsv`, which can be fed to `human_best_hits.py`:
```
python human_best_hits.py human_distances.tsv
```
