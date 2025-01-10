#!/bin/bash
set -e

while IFS="," read -r -u 3 accession tolid SRA
do
  esearch -db sra -query $SRA | esummary | xtract -pattern DocumentSummary -element Sample@acc Run@acc Experiment@acc Platform instrument_model LIBRARY_STRATEGY Summary -element Statistics@total_bases'
done 3< <(grep 'ERS\|SRS' raw_data_metadata.ls | grep -v alt)
