datasets summary genome accession $(cat zoonomia_acc.txt) --as-json-lines > zoonomia.jsonl
datasets summary genome accession $(cat zoonomia_acc.txt) --as-json-lines --reference > zoonomia_ref.jsonl
datasets summary genome accession $(cat DNAZoo_acc.txt) --as-json-lines > DNAZoo.jsonl
datasets summary genome accession $(cat DNAZoo_acc.txt) --as-json-lines --reference > DNAZoo_ref.jsonl

for g in $(cat iridian_acc.txt);do
	datasets summary genome accession $g --as-json-lines >> iridian.jsonl
	datasets summary genome accession $g --as-json-lines --reference >> iridian_ref.jsonl
done
datasets summary genome accession PRJNA545868 --as-json-lines > B10K.jsonl
datasets summary genome accession PRJNA545868 --as-json-lines --reference > B10K_ref.jsonl
datasets summary genome accession PRJNA550295 --as-json-lines > CICHLIDX.jsonl 
datasets summary genome accession PRJNA550295 --as-json-lines --reference > CICHLIDX_ref.jsonl 
