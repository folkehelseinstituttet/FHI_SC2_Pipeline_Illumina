#!/bin/bash

# Wrapper script for running "create_Gisaid_submission_files.R"
# Example run:
# ./Binaries/Gisaid_create_submission.sh

# Change the "oppsett" below. E.g. "FHI200", or "681" for Run681
declare -a array=("FHI247" "FHI249" "FHI254" "FHI256" "FHI257" "FHI259")
platform="Swift_FHI"
submitter="jonbra"

# loop through arguments
#while getopts p:o:f:m:S: flag
#do
#    case "${flag}" in
#        p) platform=${OPTARG};;
#        o) oppsett=${OPTARG};;
#        f) fasta=${OPTARG};;
#        m) metadata=${OPTARG};;
#        S) submitter=${OPTARG};;
#    esac
#done
#echo "Platform: $platform";
#echo "Oppsett: $oppsett";
#echo "Fasta filename: $fasta";
#echo "Metadata filename: $metadata";
#echo "Username submitter: $submitter";

# make tmp directory and go to it
mkdir -p ~/tmp_gisaid/Frameshift
cd ~/tmp_gisaid/

# Run R script
for oppsett in "${array[@]}"; do \
  docker run -it --rm \
  -v /mnt/N/Virologi/:/home/docker/N \
  -v $(pwd):/home/docker/Fastq \
  garcianacho/fhisc2:Illumina \
  Rscript --vanilla /home/docker/Scripts/Gisaid_create_submission_files.R -p ${platform} -o "${oppsett}" -f "${oppsett}.fasta" -m "${oppsett}.csv" -S ${submitter}; done

#head -1 ${array[1]} > date.csv
#for oppsett in "${array[@]}"; do grep -v "^covv" "${oppsett}.csv" >> date.csv; done
#for oppsett in "${array[@]}"; do cat "${oppsett}.fasta" >> date.fasta; done
