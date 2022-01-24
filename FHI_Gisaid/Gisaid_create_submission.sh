#!/bin/bash

# Wrapper script for running "create_Gisaid_submission_files.R"
# Save a copy of this file and keep it with the submission files for references.
# Example run:
# Binaries/Gisaid_create_submission.sh

# Set the name of the submission directory on N: (NB: This folder has to be created before you start)
# Example: subm_dir="2022-01-20_FHI_batch"
# subm_dir="2022-00-00_test"

# Change the "oppsett" below. 
# E.g. "FHI200" for Swift_FHI, 
# or "681" for Artic_Illumina Run681, 
# Nr134A/Nano for Artic_Nanopore, 
# or MIK172 for Swift_MIK
declare -a array=("MIK311" "MIK312" "MIK314")

# Change the platform type. Allowed values:
# Illumina NSC (FHI): Swift_FHI
# Illumina Artic: Artic_Illumina
# Illumina MIK: Swift_MIK
# Nanopore: Artic_Nanopore
platform="Swift_MIK"

# Sett inn Gisaid brukernavn:
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

# Get the latest updates
cd ~/.fhiscripts/FHI_SC2_Pipeline_Illumina
git pull origin master
docker build -t garcianacho/fhisc2:Illumina .

# make tmp directory and go to it
mkdir -p ~/tmp_gisaid/Frameshift
cd ~/tmp_gisaid/

# Clear any old submission files
rm *.{csv,fasta,xlsx}

# Run R script
for oppsett in "${array[@]}"; do \
  docker run -it --rm \
  -v /mnt/N/Virologi/:/home/docker/N \
  -v $(pwd):/home/docker/Fastq \
  garcianacho/fhisc2:Illumina \
  Rscript --vanilla /home/docker/FHI_Gisaid/Gisaid_create_submission_files.R -p ${platform} -o "${oppsett}" -f "${oppsett}.fasta" -m "${oppsett}.csv" -S ${submitter}; done

# Merge the files for each oppsett
now=`date +"%Y-%m-%d"`
declare -a array=(*.csv)
head -1 "${array[0]}" > ${now}.csv
for oppsett in "${array[@]}"; do grep -v "^covv" "${oppsett}" >> ${now}.csv; done
for oppsett in "${array[@]}"; do cat "${oppsett%csv}fasta" >> ${now}.fasta; done

# Move the submission files to N
# mv *.{csv,fasta,xlsx} /mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/4-GISAIDsubmisjon/${subm_dir}
# Move a copy of the submission script for reference
# mv ~/.fhiscripts/FHI_SC2_Pipeline_Illumina/*.sh /mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/4-GISAIDsubmisjon/${subm_dir}
