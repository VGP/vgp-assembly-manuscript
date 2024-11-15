#conda activate jq

output_file='vgp_assemblyStats.txt'
echo -e "Accession\tComplete_BUSCOs\tSingle_copy_BUSCOs\tDuplicate_BUSCOs\tFragmented_BUSCOs\tMissing_BUSCOs\tInternal_stop_codons\tTax_ID\tSubmission_Date\tNumber_Contigs\tContig_N50\tContig_N90\tNumber_Scaffolds\tScaffold_N50\tScaffold_N90\tAssembly_Size\tPercentage_masked\tSpecies_name\tSuper_class\tClass\tInfra_class\tSub_class\tSuper_order\tOrder\tSequencing_type" > "$output_file"

for gca in $(cat gca.txt);do

	busco_json=$(ls ncbi_dataset/data/${gca}/${gca}.miniprot.busco/short*json)

	complete_buscos=$(jq '.results["Complete BUSCOs"]' ${busco_json})
	single_copy_buscos=$(jq '.results["Single copy BUSCOs"]' ${busco_json})
	duplicate_buscos=$(jq '.results["Multi copy BUSCOs"]' ${busco_json})
	fragmented_buscos=$(jq '.results["Fragmented BUSCOs"]' ${busco_json})
	missing_buscos=$(jq '.results["Missing BUSCOs"]' ${busco_json})
	internal_stop_codons=$(jq '.results.internal_stop_codon_count' ${busco_json})

	tax_id=$(jq --arg acc "$gca" -r '. | select(.accession == $acc) | .organism.tax_id' genome_metadata.json)
	date=$(jq --arg acc "$gca" -r '. | select(.accession == $acc) | .assembly_info.release_date' genome_metadata.json)
	methods=$(jq --arg acc "$gca" -r '. | select(.accession == $acc) | .assembly_info.assembly_method' genome_metadata.json)
	datatypes=$(jq --arg acc "$gca" -r '. | select(.accession == $acc) | .assembly_info.sequencing_tech' genome_metadata.json)
	comments=$(jq --arg acc "$gca" -r '. | select(.accession == $acc) | .assembly_info.comments' genome_metadata.json)
	

	if echo "$methods" | grep -q '[Hh]i[Ff]i'; then
	    data='PacBio HiFi'
	elif echo "$methods" | grep -q 'IPA'; then
	    data='PacBio HiFi'
	elif echo "$methods" | grep -q '[Hh]icanu'; then
	    data='PacBio HiFi'
	elif echo "$datatypes" | grep -q '[Hh]i[Ff]i'; then
	    data='PacBio HiFi'
	elif echo "$comments" | grep -q '[Hh]i[Ff]i'; then
	    data='PacBio HiFi'
	elif echo "$datatypes" | grep -q '[Oo]xford'; then
	    data='Oxford Nanopore'
	elif echo "$datatypes" | grep -q '[Nn]anopore'; then
	    data='Oxford Nanopore'
	elif echo "$datatypes" | grep -q 'ONT'; then
	    data='Oxford Nanopore'
	elif echo "$comments" | grep -q '[Oo]xford'; then
	    data='Oxford Nanopore'
	elif echo "$comments" | grep -q '[Nn]anopore'; then
	    data='Oxford Nanopore'
	elif echo "$comments" | grep -q 'ONT'; then
	    data='Oxford Nanopore'
	elif echo "$datatypes" | grep -q '[Cc][Ll][Rr]'; then
	    data='PacBio CLR'
	elif echo "$comments" | grep -q '[Cc][Ll][Rr]'; then
	    data='PacBio CLR'
	elif echo "$methods" | grep -q '[Cc][Ll][Rr]'; then
	    data='PacBio CLR'
	elif echo "$methods" | grep -q '[Ff]alcon'; then
	    data='PacBio CLR'
	elif echo "$methods" | grep -q 'FALCON'; then
	    data='PacBio CLR'
	elif echo "$comments" | grep -q '[Ff]alcon'; then
	    data='PacBio CLR'
	elif echo "$datatypes" | grep -q '[Rr][Ss][Ii][Ii]'; then
	    data='PacBio RSII'
	else
	    data='missing'
	    echo $gca
	fi

	assembly_summary=$(ls ncbi_dataset/data/${gca}/*stats)

	scaffold_count=$(cat "$assembly_summary" | grep "# scaffolds:" | awk '{print $3}')
	total_scaffold_length=$(cat "$assembly_summary" | grep "Total scaffold length:" | awk '{print $4}')
	scaffold_n50=$(cat "$assembly_summary" | grep "Scaffold N50:" | head -n 1 | awk '{print $3}')
	scaffold_n90=$(cat "$assembly_summary" | grep "Scaffold N90:" | awk '{print $3}')
	contig_count=$(cat "$assembly_summary" | grep "# contigs:" | awk '{print $3}')
	total_contig_length=$(cat "$assembly_summary" | grep "Total contig length:" | awk '{print $4}')
	contig_n50=$(cat "$assembly_summary" | grep "Contig N50:" | head -n 1 | awk '{print $3}')
	contig_n90=$(cat "$assembly_summary" | grep "Contig N90:" | awk '{print $3}')
	soft_masked_bases=$(cat "$assembly_summary" | grep "# soft-masked bases:" | awk '{print $4}')
	
	# Calculating total length
	total_length=$((total_contig_length))
	
	# Calculating percentage of genome soft-masked
	soft_masked_percentage=$(awk "BEGIN {printf \"%.2f\", ($soft_masked_bases / $total_length) * 100}")

	curl -s https://www.ebi.ac.uk/ena/browser/api/xml/${tax_id} -o response.xml
	species_name=$(xmllint --xpath 'string(//taxon/@scientificName)' response.xml)

	# Extract the entire lineage
	superclass=$(grep 'rank="superclass"' response.xml | awk -F 'scientificName="' '{print $2}' | awk -F '"' '{print $1}')
	class=$(grep 'rank="class"' response.xml | awk -F 'scientificName="' '{print $2}' | awk -F '"' '{print $1}')
	infraclass=$(grep 'rank="infraclass"' response.xml | awk -F 'scientificName="' '{print $2}' | awk -F '"' '{print $1}')
	subclass=$(grep 'rank="subclass"' response.xml | awk -F 'scientificName="' '{print $2}' | awk -F '"' '{print $1}')
	superorder=$(grep 'rank="superorder"' response.xml | awk -F 'scientificName="' '{print $2}' | awk -F '"' '{print $1}')
	order=$(grep 'rank="order"' response.xml | awk -F 'scientificName="' '{print $2}' | awk -F '"' '{print $1}')


	echo -e "${gca}\t${complete_buscos}\t${single_copy_buscos}\t${duplicate_buscos}\t${fragmented_buscos}\t${missing_buscos}\t${internal_stop_codons}\t${tax_id}\t${date}\t${contig_count}\t${contig_n50}\t${contig_n90}\t${scaffold_count}\t${scaffold_n50}\t${scaffold_n90}\t${total_scaffold_length}\t${soft_masked_percentage}\t${species_name}\t${superclass}\t${class}\t${infraclass}\t${subclass}\t${superorder}\t${order}\t${data}" >> "$output_file"
done

