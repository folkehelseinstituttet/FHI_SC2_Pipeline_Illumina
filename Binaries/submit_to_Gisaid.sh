#!/bin/bash
# Run this script after someone has controlled and approved a Gisaid create_Gisaid_submission_files

#cd ~/tmp_gisaid/

# Declare an array of the "oppsett" you want to submit. E.g.:
#declare -a array=("FHI200" "FHI201" "FHI202")

# Create directories on N:
#N_path="/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/4-GISAIDsubmisjon/"

#for oppsett in "${array[@]}"; do sudo mkdir ${N_path}$(date +"%Y-%m-%d")_${oppsett}; done

# Move files to N:
#for oppsett in "${array[@]}"; do sudo cp *${oppsett}.* ${N_path}$(date +"%Y-%m-%d")_${oppsett}; done

# Her kan jeg fortsette i bash og flytte på filer til N?, og kanskje ta gisaid-submission?
# Alle filene jeg trenger blir liggende i tmp_gisaid-mappenavn
# Fra her kan jeg kjøpre cli2-scriptet
# cli2 upload --token /home/jonr/Downloads/jonbra_gisaid_token \
#  --metadata FHI203.csv \
#  --fasta FHI203.fasta \
#  --database EpiCoV \
#  --frameshift catch_none \
#  --failed FHI203_failed_upload.csv \
#  --log FHI203_submisjon.log

#sudo mv FHI203_failed_upload.csv FHI203_submisjon.log
