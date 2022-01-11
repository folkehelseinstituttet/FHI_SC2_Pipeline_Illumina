#!/bin/bash

# Wrapper script for running "create_Gisaid_submission_files.R"
# Save a copy of this file and keep it with the submission files for references.
# Example run:
# ./Gisaid_create_submission.sh

# Change the "oppsett" below. E.g. "FHI200" (FHI NSC), or "681" (Artic Illumina Run681), Nr134A/Nano (Nanopore), MIK172 (MIK NSC)
declare -a array=("FHI271a" "FHI274" "FHI275" "FHI276" "FHI277")

# Illumina NSC (FHI): Swift_FHI
# Illumina Artic: Illumina_Artic
# Illumina MIK: Swift_MIK
# Nanopore: Artic_Nanopore
platform="Swift_FHI"

# Sett inn Gisaid brukernavn:
submitter="hildenf"

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

now=`date +"%Y-%m-%d"`
declare -a array=(*.csv)
head -1 "${array[0]}" > ${now}.csv
for oppsett in "${array[@]}"; do grep -v "^covv" "${oppsett}" >> ${now}.csv; done
for oppsett in "${array[@]}"; do cat "${oppsett%csv}fasta" >> ${now}.fasta; done
#head -1 ${array[1]} > date.csv
#for oppsett in "${array[@]}"; do grep -v "^covv" "${oppsett}.csv" >> date.csv; done
#for oppsett in "${array[@]}"; do cat "${oppsett}.fasta" >> date.fasta; done
