#!/bin/bash
basedir=/home/docker/Fastq
runname=${1%/}
runname=${runname##/*/}

cat /home/docker/mountingpoint/*.fa*  > /home/docker/Fastq/${runname}.fa

#RefLength=29903
seqkit fx2tab -C A,C,T,G,N ${runname}.fa | awk '{print $1,($3+$4+$5+$6)/29903*100}' > data_file.csv
sed 's/\ /,/g' data_file.csv > data_file2.csv
echo "name2,coverage_breath" > header_file.csv
cat header_file.csv data_file2.csv > coverage.csv


source activate pangolin
pangolin --update
pangolin ./${runname}.fa -t 8 --outfile ./${runname}_pangolin_out.csv
conda deactivate

nextclade --input-fasta /home/docker/Fastq/${runname}.fa  --output-csv /home/docker/Fastq/${runname}_Nextclade.results.csv
nextalign  --sequences=${basedir}/${runname}_summaries/fasta/${runname}.fa --reference=/home/docker/CommonFiles/reference_nc.fasta \
 --genemap=/home/docker/CommonFiles/genemap.gff --genes=E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S --output-dir=${basedir} --output-basename=${runname}
Rscript /home/docker/Scripts/InsertionAnalysis.R

nextclade_output_converter.py ${runname}_Nextclade.results.csv >> ${runname}_Nextclade.results2.csv
awk -F ',' '{print $1 "," $2 "," $4}' ${runname}_pangolin_out.csv > pangolin_out.csv
awk -F ';' '{print $1 "," $2}' ${runname}_Nextclade.results.csv > nextclade_out2.csv

cat nextclade_out2.csv | sed "s/, /\//g" > nextclade_out3.csv && mv nextclade_out3.csv nextclade_out2.csv #ny fra 22.06.21 Kamilla&Nacho

(head -n 1 pangolin_out.csv && tail -n +2 pangolin_out.csv | sort) > pangolin_out_sorted.csv
(head -n 1 nextclade_out2.csv && tail -n +2 nextclade_out2.csv | sort) > nextclade.out2_sorted.csv
(head -n 1 ${runname}_Nextclade.results2.csv && tail -n +2 ${runname}_Nextclade.results2.csv | sort) > ${runname}_Nextclade.results2_sorted.csv
(head -n 1 coverage.csv && tail -n +2 coverage.csv | sort) > coverage2.csv

paste -d, ${runname}_Nextclade.results2_sorted.csv pangolin_out_sorted.csv > NextcladeAndPangolin.out.csv
paste -d, NextcladeAndPangolin.out.csv nextclade.out2_sorted.csv > NextcladeAndPangolin.out2.csv
paste -d, NextcladeAndPangolin.out2.csv coverage2.csv > NextcladeAndPangolinAndCoverage.csv

sed 's/,/\t/g' NextcladeAndPangolinAndCoverage.csv | sed 's/ORF10/ORF10\t/g' > ${runname}_NextcladeAndPangolin.csv


rm *results*
rm *out*
rm *file*
rm *age*
rm ./${runname}.fa
rm /home/docker/Fastq/*.fasta
rm /home/docker/Fastq/*.insertions.csv
rm /home/docker/Fastq/*.errors.csv

echo "Pangolin og Nextclade ferdig"
