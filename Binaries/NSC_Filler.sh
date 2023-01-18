#!/bin/bash

source activate nextclade
nextclade dataset get --name 'sars-cov-2' --output-dir '/home/docker/nc_sars-cov-2'
nextclade --input-fasta ivar_masked_consensus_Nremoved_v*.fa --input-dataset /home/docker/nc_sars-cov-2 --output-csv Nextclade.new.results.csv
conda deactivate

Rscript /home/docker/Scripts/LongPangoParserNSC.R

rm Nextclade.new.results.csv
rm *.gene.*.fasta
rm *.aligned.fasta
rm *.errors.csv
rm *.insertions.csv
