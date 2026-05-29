for f in ncbi_dataset/data/*/*fna.gz;do
	seqkit fx2tab -n -l $f > ${f%.fna.gz}_scaffolds.tsv
	seqkit seq -w 0 ${f} | \
		awk '/^>/{if(seq && hdr){seq=toupper(seq); m=split(seq,p,/N+/); for(i=1;i<=m;i++) if(length(p[i])>0) printf("%s_contig%06d\t%d\n",hdr,i,length(p[i]))}; hdr=substr($0,2); seq=""; next} {seq=seq $0} END{if(seq && hdr){seq=toupper(seq); m=split(seq,p,/N+/); for(i=1;i<=m;i++) if(length(p[i])>0) printf("%s_contig%06d\t%d\n",hdr,i,length(p[i]))}}' > ${f%.fna.gz}_contigs.tsv
done
