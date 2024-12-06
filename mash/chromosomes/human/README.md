# Create a lookup table for human chromosome names:
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/695/395/GCA_016695395.2_mHomSap3.mat/GCA_016695395.2_mHomSap3.mat_assembly_report.txt
grep "^chr" GCA_016695395.2_mHomSap3.mat_assembly_report.txt | cut -f1,5 > chr.maternal.lookup.tsv
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/700/455/GCA_016700455.2_mHomSap3.pat/GCA_016700455.2_mHomSap3.pat_assembly_report.txt
grep "^chr" GCA_016700455.2_mHomSap3.pat_assembly_report.txt | cut -f1,5 >> chr.paternal.lookup.tsv
cat chr.maternal.lookup.tsv chr.paternal.lookup.tsv | grep -f <(gfastats -ss complete.fasta | cut -f1) > chr.combined.lookup.tsv
```

# Collect sequences (including x,y chrs):

```
datasets download genome accession GCA_016695395.2 --filename GCA_016695395.2.zip
unzip -o GCA_016695395.2.zip -d GCA_016695395.2
datasets download genome accession GCA_016700455.2 --filename GCA_016700455.2.zip
unzip -o GCA_016700455.2.zip -d GCA_016700455.2
cat GCA_016695395.2/*/*/*/*.fna GCA_016700455.2/*/*/*/*.fna > combined.fasta
gfastats -i <(cut -f2 chr.combined.lookup.tsv) -o complete.fasta
```

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
You can further assign human chromosome numbers using `assign_chr.py`:
```
assign_chr.py human_outliers.csv chr.combined.lookup.tsv 
```