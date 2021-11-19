#!/bin/bash

# Wrapper script for running "create_Gisaid_submission_files.R"

# loop through arguments
while getopts p:o:f:m:S: flag
do
    case "${flag}" in
        p) platform=${OPTARG};;
        o) oppsett=${OPTARG};;
        f) fasta=${OPTARG};;
        m) metadata=${OPTARG};;
        S) submitter=${OPTARG};;
    esac
done
echo "Platform: $platform";
echo "Oppsett: $oppsett";
echo "Fasta filename: $fasta";
echo "Metadata filename: $metadata";
echo "Username submitter: $submitter";

# make tmp directory and go to it
mkdir -p ~/tmp_gisaid/Frameshift
cd ~/tmp_gisaid/

# Run R script
docker run -it --rm \
  -v /mnt/N/Virologi/:/home/docker/N \
  -v $(pwd):/home/docker/Fastq \
  garcianacho/fhisc2:Illumina \
  Rscript --vanilla /home/docker/Scripts/create_Gisaid_submission_files.R -p $platform -o $oppsett -f $fasta -m $metadata -S $submitter

# Her kan jeg fortsette i bash og flytte på filer til N?, og kanskje ta gisaid-submission?