(
  echo -e "accession\ttaxid\tspecies_name\tassembly_level\tmethod\tsequencing_technology\tsubmission_date\tnumber_of_contigs\tnumber_of_scaffolds\tcontig_n50\tscaffold_n50"
  jq -r '
    [
      .accession,
      .organism.tax_id,
      .organism.organism_name,
      .assembly_info.assembly_level,
      .assembly_info.assembly_method,
      .assembly_info.sequencing_tech,
      .assembly_info.biosample.submission_date,
      .assembly_stats.number_of_contigs,
      .assembly_stats.number_of_scaffolds,
      .assembly_stats.contig_n50,
      .assembly_stats.scaffold_n50
    ] | @tsv
  ' DNAZoo_ref.jsonl
) > DNAZoo_ref.metadata.tsv

(
  echo -e "accession\ttaxid\tspecies_name\tassembly_level\tmethod\tsequencing_technology\tsubmission_date\tnumber_of_contigs\tnumber_of_scaffolds\tcontig_n50\tscaffold_n50"
  jq -r '
    [
      .accession,
      .organism.tax_id,
      .organism.organism_name,
      .assembly_info.assembly_level,
      .assembly_info.assembly_method,
      .assembly_info.sequencing_tech,
      .assembly_info.biosample.submission_date,
      .assembly_stats.number_of_contigs,
      .assembly_stats.number_of_scaffolds,
      .assembly_stats.contig_n50,
      .assembly_stats.scaffold_n50
    ] | @tsv
  ' zoonomia_ref.jsonl
) > zoonomia_ref.metadata.tsv

(
  echo -e "accession\ttaxid\tspecies_name\tassembly_level\tmethod\tsequencing_technology\tsubmission_date\tnumber_of_contigs\tnumber_of_scaffolds\tcontig_n50\tscaffold_n50"
  jq -r '
    [
      .accession,
      .organism.tax_id,
      .organism.organism_name,
      .assembly_info.assembly_level,
      .assembly_info.assembly_method,
      .assembly_info.sequencing_tech,
      .assembly_info.biosample.submission_date,
      .assembly_stats.number_of_contigs,
      .assembly_stats.number_of_scaffolds,
      .assembly_stats.contig_n50,
      .assembly_stats.scaffold_n50
    ] | @tsv
  ' B10K_ref.jsonl
) > B10K_ref.metadata.tsv

(
  echo -e "accession\ttaxid\tspecies_name\tassembly_level\tmethod\tsequencing_technology\tsubmission_date\tnumber_of_contigs\tnumber_of_scaffolds\tcontig_n50\tscaffold_n50"
  jq -r '
    [
      .accession,
      .organism.tax_id,
      .organism.organism_name,
      .assembly_info.assembly_level,
      .assembly_info.assembly_method,
      .assembly_info.sequencing_tech,
      .assembly_info.biosample.submission_date,
      .assembly_stats.number_of_contigs,
      .assembly_stats.number_of_scaffolds,
      .assembly_stats.contig_n50,
      .assembly_stats.scaffold_n50
    ] | @tsv
  ' CICHLIDX_ref.jsonl
) > CICHLIDX_ref.metadata.tsv

(
  echo -e "accession\ttaxid\tspecies_name\tassembly_level\tmethod\tsequencing_technology\tsubmission_date\tnumber_of_contigs\tnumber_of_scaffolds\tcontig_n50\tscaffold_n50"
  jq -r '
    [
      .accession,
      .organism.tax_id,
      .organism.organism_name,
      .assembly_info.assembly_level,
      .assembly_info.assembly_method,
      .assembly_info.sequencing_tech,
      .assembly_info.biosample.submission_date,
      .assembly_stats.number_of_contigs,
      .assembly_stats.number_of_scaffolds,
      .assembly_stats.contig_n50,
      .assembly_stats.scaffold_n50
    ] | @tsv
  ' ebp_other_ref.jsonl
) > ebp_other_ref.metadata.tsv

(
  echo -e "accession\ttaxid\tspecies_name\tassembly_level\tmethod\tsequencing_technology\tsubmission_date\tnumber_of_contigs\tnumber_of_scaffolds\tcontig_n50\tscaffold_n50"
  jq -r '
    [
      .accession,
      .organism.tax_id,
      .organism.organism_name,
      .assembly_info.assembly_level,
      .assembly_info.assembly_method,
      .assembly_info.sequencing_tech,
      .assembly_info.biosample.submission_date,
      .assembly_stats.number_of_contigs,
      .assembly_stats.number_of_scaffolds,
      .assembly_stats.contig_n50,
      .assembly_stats.scaffold_n50
    ] | @tsv
  ' iridian_ref.jsonl
) > iridian_ref.metadata.tsv

(
  echo -e "accession\ttaxid\tspecies_name\tassembly_level\tmethod\tsequencing_technology\tsubmission_date\tnumber_of_contigs\tnumber_of_scaffolds\tcontig_n50\tscaffold_n50"
  jq -r '
    [
      .accession,
      .organism.tax_id,
      .organism.organism_name,
      .assembly_info.assembly_level,
      .assembly_info.assembly_method,
      .assembly_info.sequencing_tech,
      .assembly_info.biosample.submission_date,
      .assembly_stats.number_of_contigs,
      .assembly_stats.number_of_scaffolds,
      .assembly_stats.contig_n50,
      .assembly_stats.scaffold_n50
    ] | @tsv
  ' vertebrata_ref.jsonl
) > vertebrata_ref.metadata.tsv

