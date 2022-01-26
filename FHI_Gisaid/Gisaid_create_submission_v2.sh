#!/bin/bash

# Wrapper script for running "Gisaid_create_submission_files.R"
# Example run (while in .fhiscripts folder):
# ./Gisaid_create_submission_v2.sh -f "2022-01-26-test/Gisaid_sample_sheet.xlsx"
# "2022-01-26-test" refers to the submission directory on the N drive

# loop through arguments
while getopts ":f:" flag; do
    case "${flag}" in
        f) 
          folder=${OPTARG}
          ;;
    esac
done

# Get the latest script updates from Github
cd ~/.fhiscripts/FHI_SC2_Pipeline_Illumina
git pull origin master
docker build -t garcianacho/fhisc2:Illumina .

# make a tmp directory and go to it
mkdir -p ~/tmp_gisaid/Frameshift
cd ~/tmp_gisaid/

# Clear any old submission files
rm *.{csv,fasta,log}
rm FrameShift*.xlsx

# Run the submission script
docker run -it --rm \
  -v /mnt/N/Virologi/:/home/docker/N \
  -v $(pwd):/home/docker/Fastq \
  garcianacho/fhisc2:Illumina Rscript --vanilla /home/docker/FHI_Gisaid/Gisaid_create_submission_files_v2.R ${folder}
