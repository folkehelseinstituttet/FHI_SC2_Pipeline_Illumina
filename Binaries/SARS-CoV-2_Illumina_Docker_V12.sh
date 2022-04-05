#!/bin/bash
#### NY VERSJON V2 ######   # Byttet primer-fil fra V1 til V3 
                            # fordi filnavn er byttet til å starte på Artic istedenfor Virus er for-loopene oppdatert

#### NY VERSJON V3 ######   #bcftools mpileup d 100000
                            #erstatttet variant calling og consensus med ivar
#### NY VERSJON V4 ######   #v4b: Q>30
                            #v4c La til B flagget. Anbefalt i iVar manual. (Disable base alignment quality (BAQ) computation. See BAQ below.) 
                            #v4c.t0  majority rules
#### NY VERSJON V5 ######   #Lagt til automatisk kopiering av resultatfiler til felles. År og DirFelles må byttes ut avhengig av årstall og linux-maskin. Må skrive inn passord til linux-maskinen i terminalen for å starte kopieringen. 

#### NY VERSJON V6 ######   #Lagt til fjerning av N'er i starten og slutten av fasta-fil (consensus)

#### NY VERSJON V7 ######   #ivar oppdatert fra versjon 1.0 til versjon 1.3
                            #Lagt til pangolin og nextclade, inkl. samle alle konsensus  i en fasta
                            #Endret filbane summaries kopieres til /mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina/

#### NY VERSJON V8 ######   #Rettet opp forskyvning av SeqID på pangolin-resultater

#### NY VERSJON V9 ######   #Lagt til variabel for linux-maskin og årstall
							#Oppdatert hvor på N: summaries og fasta-fil for hele runnet kopieres til
							#Lagt til kolonnen for pangolin versjon i output
							#Fikset trim-rapport til å passe med ivar versjon 1.3
							#Samlet *NextcladeAndPangolin.csv og *summaries.csv til *summaries_and_Pangolin.csv
#VERSJON V9.1 liten justering#Lagt til V4 av bed-fil med primeroversikt

#### NY VERSJON V10 ######  #Sletter overflødige csv fra N: etter kopiering av summaries-mappen
							#Lagt til R-skript som lager fasta med kun spike-sekvensene


# henter ut året for 4 dager siden  (ved årsskiftet vil et run kunne tilhøre året før)
# Docker version of script

year=$([ "$OSTYPE" = linux-gnu ] && date --date="4 days ago" +"%Y" || date -v-4d +"%Y")	                       

basedir=$(pwd)
runname=$(ls *.xlsx)
runname=${runname%.xlsx}

cd /home/docker/Fastq

if [ ${1} == "ArticV3" ]; then
cp /home/docker/CommonFiles/artic-ncov2019/primer_schemes/nCoV-2019/V3.2/nCoV-2019.primer.bed /home/docker/Fastq/primers.bed
STR="/home/docker/CommonFiles/artic-ncov2019/primer_schemes/nCoV-2019/V3.2/nCoV-2019.primer.bed /home/docker/Fastq/primers.bed"
primer_version="${STR:65:4}"

fi

if [ ${1} == "ArticV4" ]; then
cp /home/docker/CommonFiles/artic-ncov2019/primer_schemes/nCoV-2019/V4.1/SARS-CoV-2.primer.bed /home/docker/Fastq/primers.bed
STR="/home/docker/CommonFiles/artic-ncov2019/primer_schemes/nCoV-2019/V4.1/nCoV-2019.primer.bed /home/docker/Fastq/primers.bed"
primer_version="${STR:65:4}"

fi


if [[ ${1} ]]; then  
if [ ${1} != "ArticV3" ]  && [ ${1} != "ArticV4" ]; then
       echo "
  _  _ ___ _ 
 | \| | _ ) |
 | .  | _ \_|
 |_|\_|___(_)

            NÅ SKREV DU FEIL 
            Du skrev ${1} som input istedenfor Artic4 eller Artic3 
#############################################################################"
    else  
    echo "       
                    ################################
                  ##                                ##
                ##                                   ##
               ##   Du bruker nå primersett           ##
               ##  ${1}             ##
               ##                                    ##
                ##                                  ##
                  ##                               ##  
                   ################################"
fi
else
   echo "
  _  _ ___ _ 
 | \| | _ ) |
 | .  | _ \_|
 |_|\_|___(_)
 
            Du har glemt å skrive ArticV3, ArticV4 eller Midnight 
#############################################################################  


"
fi


cp /home/docker/CommonFiles/nCoV-2019.reference.fasta /home/docker/Fastq
#cp /home/docker/CommonFiles/artic-ncov2019/primer_schemes/nCoV-2019/V3.2/nCoV-2019.bed /home/docker/Fastq/primers.bed
scriptdir=/home/docker/Scripts
script_name1=`basename $0`

###### DATABASER/REFERANSESEKVENSER ########

CoronaRef=/home/docker/Fastq/nCoV-2019.reference.fasta
PrimerBed=/home/docker/Fastq/primers.bed

########## FYLL INN FOR AGENS ###################
#Skriv inn agens-navn (må være skrevet likt som i navnet på fasta-fil som inneholder databasen/referansesekvensene
Agens=Corona					#ingen mellomrom etter =
#presisere stringency for mapping, 1-100
String2=85 

######## DEL 1 Mapping #### START ######

for dir in $(ls -d */)
do
    cd ${dir}

	R1=$(ls *_R1*.fastq.gz)
    R2=$(ls *_R2*.fastq.gz)

gzip -d ${R1} 
gzip -d ${R2} 

    newR1=$(ls *_R1*.fastq)
    newR2=$(ls *_R2*.fastq)
#align
    bestF2="${newR1%%_*L001*}_" # brukes til navnsetting av outputfil 
    tanoti -r ${CoronaRef} -i ${newR1} ${newR2} -o ${bestF2}tanoti.sam -p 1 -m ${String2}
    bestF3=$(ls *_tanoti.sam)
    samtools view -bS ${bestF3} | samtools sort -o ${bestF3%.sam}_sorted.bam
    samtools index ${bestF3%.sam}_sorted.bam
   
    cd "${basedir}"
done

echo "HEY HEY HEY, What's that sound?" 
echo "Mapping done!"

######## DEL 1 Mapping #### SLUTT ######

######## DEL 2 Trimming #### START ######
#source activate ivar13
#we can now install ivar13 in corona environment#

cd "${basedir}"

for dir in $(ls -d */)
do
    cd ${dir}
    bestF4=$(ls *_sorted.bam)
#trim
    ivar trim -i ${bestF4} -b ${PrimerBed} -p ${bestF4%_sorted.bam}.trimmed -m 50 -q 15 -s 4 > ${bestF4%_sorted.bam}.TrimReport.txt #brukt forslag fra manual
    
    samtools view -bS ${bestF4%_sorted.bam}.trimmed.bam | samtools sort -o ${bestF4%_sorted.bam}.trimmed.sorted.bam
    samtools index ${bestF4%_sorted.bam}.trimmed.sorted.bam
    cd "${basedir}"
done
#conda deactivate

echo "#
"
echo "Read ferdig trimmet" 
echo "

"
echo "###################"
######## DEL 2 Trimming #### SLUTT ######

######## DEL 3 VariantCalling og Consensus #### START ######

cd "${basedir}"

for dir in $(ls -d Art*/)
do

cd ${dir}

# Lage konsensus for Main-genotype
    bestF5=$(ls *trimmed.sorted.bam)

	samtools sort -n ${bestF5} > ${bestF5%.bam}.byQuery.bam 
#tester samtools vs bcftools. Versjon 1.9
	samtools mpileup -f ${CoronaRef} ${bestF5%} -d 100000 -Q 30 -B| ivar consensus -t 0 -m 10 -n N -p cons.fa 
    #bcftools index calls.vcf.gz
    #bedtools genomecov -bga -ibam ${bestF5}| grep -w '0$\|1$\|2$\|3$' > regionswithlessthan4coverage.bed   
    #ivar consensus -m regionswithlessthan4coverage.bed -f ${CoronaRef} calls.vcf.gz -o cons.fa
    seqkit replace -p "(.+)" -r ${bestF5%%_*} cons.fa | sed -r '2s/^N{1,}//g' | sed -r '$ s/N{1,}$//g' > ${bestF5%%_*}_consensus.fa

#sletter filer som ikke trengs videre: 
	rm *cons.fa 
	#rm *calls*.vcf.gz
	#rm *calls*.vcf.gz.csi 
	#rm *regionswith*coverage.bed 
	rm *sorted.byQuery.bam 

cd "${basedir}"

done

echo "
ivar1.3 consensus done
"
######## DEL 3 VariantCalling og Consensus #### SLUTT ######

######## DEL 4 CoveragePlot og Statistikk #### START ######
cd "${basedir}"

for dir in $(ls -d Art*/)
do
    cd ${dir}
	
# Coverage plot og statistik 
	#etter trimming
	bestF5=$(ls *trimmed.sorted.bam)
	weeSAMv1.6 --bam ${bestF5} --out ${bestF5%.bam}_stats.txt 
	bestF6=$(ls *_sorted.bam)
	weeSAMv1.6 --bam ${bestF6} --out ${bestF6%.bam}_stats.txt 
	
cd "${basedir}"

done

######## DEL 4 CoveragePlot og Statistikk #### SLUTT ######

######## DEL 5 Identifisere parametere og lage summary for hver prøve  #### START ######
cd "${basedir}"

# Går inn i hver mappe og identifiserer ulike parametere og legger det inn i en csv fil 
for dir in $(ls -d Art*/)
do
    cd ${dir}
#identify & log
    
    newR1=$(ls *_R1*.fastq)
    newR2=$(ls *_R1*.fastq)
    newR4=$(ls *_tanoti.sam)
      
   	bestF2="${newR1%%_*L001*}_"
    bestF5=$(ls *trimmed.sorted.bam)
    TrimReport=$(ls *TrimReport.txt)    

    readsb4=$(echo $(wc -l ${newR1}|awk '{print $1'}))
    readsb42=$(echo "${readsb4}/2"|bc)  #delt på 2 istedenfor 4 for å få reads for R1 og R2 (fastq fil inneholder 4 linjer per read)
    
	primerTrimed=$(grep "primers from" ${TrimReport} | awk '{print $4}')
	readstoshortlenght=$(grep "trimmed below" ${TrimReport} | awk '{print $1}' | sed 's/\%//g')
    readsoutsideprimer=$(grep "outside" ${TrimReport} | awk '{print $1}' | sed 's/\%//g')
	percReadstrim=$(echo "${readstoshortlenght}+${readsoutsideprimer}"|bc)
          
    mappedafter=$(sort -t$'\t' -k3 -nr *trimmed.*_stats.txt | grep -m1 "" | cut -f3)
    mappedb4=$(sort -t$'\t' -k3 -nr *sorted_stats.txt | grep -m1 "" | cut -f3)
   
    mapreadsper=$(echo "scale=2 ; ($mappedb4 / $readsb42) *100" | bc) #før trimming

    #wee1114=$(sort -t$'\t' -k3 -nr *_stats.txt | grep -m1 "" | cut -f5) #ny versjon av weeSAM gir 100% dekning uansett
    wee1115=$(sort -t$'\t' -k3 -nr *_stats.txt | grep -m1 "" | cut -f8)
       
	# bedtools genomecov -ibam ${bestF5} -bga > ${bestF5%.bam}_aln.bam 
    # LengthBelowDepth1=$(awk '$4 <1' *_aln.bam | awk '{a=$3-$2;print $0,a+1;}' | awk '{print $5}' | paste -sd+ | bc) #dette blir feil
    # LengthBelowDepth10=$(awk '$4 <10' *_aln.bam | awk '{a=$3-$2;print $0,a+1;}' | awk '{print $5}' | paste -sd+ | bc) #dette blir feil
    # LengthBelowDepth30=$(awk '$4 <30' *_aln.bam | awk '{a=$3-$2;print $0,a+1;}' | awk '{print $5}' | paste -sd+ | bc) #dette blir feil
    
    #RefLength=29903
    RefLength=$(awk 'FNR == 2 {print $2}' *trimmed.*_stats.txt)
    coverage_breath=$(seqkit fx2tab -C A,C,T,G,N *.fa | awk '{print ($3+$4+$5+$6)/29903*100}') #NB illumina og nanopore er ikke like her

    wee1114=$(echo "scale=5;(($RefLength-$LengthBelowDepth1)/$RefLength)*100" |bc) 
    PercCovAboveDepth9=$(echo "scale=5;(($RefLength-$LengthBelowDepth10)/$RefLength)*100" |bc)    
    PercCovAboveDepth29=$(echo "scale=5;(($RefLength-$LengthBelowDepth30)/$RefLength)*100" |bc)
  #test
#write bit
echo "Parameters, ${dir%/}" >> ${dir%/}_summary.csv
echo "Total_number_of_reads_before_mapping_and_trim:, ${readsb42}"  >> ${dir%/}_summary.csv
echo "Total_number_of_mapped_reads_befor_trim:, ${mappedb4}" >> ${dir%/}_summary.csv
echo "Percent_mapped_reads_before_trim:, ${mapreadsper}" >> ${dir%/}_summary.csv
echo "Percent_reads_trimmed_removed:, ${percReadstrim}" >> ${dir%/}_summary.csv
echo "Total_mapped_${Agens}_reads_after_trim:, ${mappedafter}" >> ${dir%/}_summary.csv
echo "Percent_covered:, n.a. " >> ${dir%/}_summary.csv
echo "Average_depth:, ${wee1115}" >> ${dir%/}_summary.csv
echo "Percent_covered_above_depth=9:, ${coverage_breath}" >> ${dir%/}_summary.csv
echo "Percent_covered_above_depth=29:, n.a. " >> ${dir%/}_summary.csv
echo "Script_name:, ${script_name1}/${primer_version}" >> ${dir%/}_summary.csv

    cd "${basedir}"
done

######## DEL 5 Identifisere parametere og lage summary for hver prøve #### SLUTT ######

######## DEL 6 Sammenfatte resultater #### START ######
cd "${basedir}"
runname=$(ls *.xlsx)
runname=${runname%.xlsx}

mkdir "${runname}_summaries"
mkdir "${runname}_summaries/fasta"
mkdir "${runname}_summaries/bam"
#mkdir "./${runname}_IGV_bam_filer"
#mkdir "./${runname}_summaries/Konsensus-sekvenser"

for dir in $(ls -d Art*/)
do

	cp ${dir}/*_summary.csv "./${runname}_summaries/"
#	cp ${dir}/*.html "./${runname}_summaries/"
#	cp ${dir}/*analysis.pdf "./${runname}_summaries/"   

#	cp ${dir}/*_consensus.fa "./${runname}_summaries/Konsensus-sekvenser/"  
	cp ${dir}/*_consensus.fa "./${runname}_summaries/fasta"
	
	cp ${dir}/*trimmed.sorted.bam "./${runname}_summaries/bam"
	
done

#lager en fil .tmp for hver prøve hvor alle verdiene legges inn i (uten overskriftene) 
	cd "./${runname}_summaries"

	for f in $(ls *.csv) 
	do
		sed 's/\./,/g' $f | awk 'BEGIN {OFS=","} {print $2}' > $f-5.tmp
    done

echo "Parameters:" >> parameters                            # Lager en fil parameteres hvor alle oversikriftene legges 
echo "Total number of reads before mapping and trim:"  >> parameters
echo "Total number of mapped reads before trim:" >> parameters
echo "Percent mapped before trimming:" >> parameters
echo "Percent reads removed with trimming:" >> parameters
echo "Total of mapped ${Agens} reads after trim:" >> parameters
echo "Percent covered:" >> parameters
echo "Average depth:" >> parameters
echo "Percent covered above depth=9:" >> parameters
echo "Percent covered above depth=29:" >> parameters

echo "Script name:" >> parameters

paste parameters *.tmp >> ${runname}_summaries.csv    # verdiene og overskriftene limes inn i en og samme fil

find . -type f -name "*.tmp" -exec rm -f {} \;
find . -type f -name "parameters" -exec rm -f {} \; # sletter de midlertidige filene
rm *summary.csv

cd "${basedir}/${runname}_summaries/fasta"
cat *.fa  > ${runname}.fa 
cp ${runname}.fa ${basedir}/${runname}_summaries/${runname}.fa #Fasta file copy to main folder

cd "${basedir}"

######## DEL 6 Sammenfatte resultater #### SLUTT ######

####### Index bam-filer i summaries-mappen ##### START #####
cd "./${runname}_summaries/bam"
for f in $(ls *sorted.bam) 
do
samtools index $f
done
cd "${basedir}"
####### Index filer i IGV-mappen ##### SLUTT #####

####### Pangolin og Nextclade  ##### START #####
source activate pangolin #Added 10Nov2021 for UShER
pangolin --update
pangolin ${basedir}/${runname}_summaries/fasta/${runname}.fa -t 8 --outfile ${basedir}/${runname}_summaries/${runname}_pangolin_out.csv 
conda deactivate #Added 10Nov2021 for UShER

nextclade --input-fasta ${basedir}/${runname}_summaries/fasta/${runname}.fa --output-csv ${basedir}/${runname}_summaries/${runname}_Nextclade.results.csv
nextalign  --sequences=${basedir}/${runname}_summaries/fasta/${runname}.fa --reference=/home/docker/CommonFiles/reference_nc.fasta \
 --genemap=/home/docker/CommonFiles/genemap.gff --genes=E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S --output-dir=${basedir} --output-basename=${runname}
Rscript /home/docker/Scripts/SpikeMissing.R
mv ${basedir}/${runname}_summaries/fasta/${runname}_Nextclade.results.csv ${basedir}/${runname}_summaries/
mv /home/docker/Fastq/MissingAA.Spike.xlsx ${basedir}/${runname}_summaries/${runname}_MissingAA.Spike.xlsx
cd "${basedir}/${runname}_summaries/"

nextclade_output_converter.py ${runname}_Nextclade.results.csv >> ${runname}_Nextclade.results2.csv

awk -F ',' '{print $1 "," $2 "," $4}' ${runname}_pangolin_out.csv > pangolin_out.csv
awk -F ';' '{print $1 "," $2}' ${runname}_Nextclade.results.csv > nextclade_out2.csv

cat nextclade_out2.csv | sed "s/, /\//g" > nextclade_out3.csv && mv nextclade_out3.csv nextclade_out2.csv #ny fra 22.06.21 Kamilla&Nacho

#(head -n 2 pangolin_out.csv && tail -n +3 pangolin_out.csv | sort) > pangolin_out_sorted.csv
(head -n 1 pangolin_out.csv && tail -n +2 pangolin_out.csv | sort) > pangolin_out_sorted.csv
#(head -n 2 nextclade_out2.csv && tail -n +3 nextclade_out2.csv | sort) > nextclade.out2_sorted.csv
(head -n 1 nextclade_out2.csv && tail -n +2 nextclade_out2.csv | sort) > nextclade.out2_sorted.csv
(head -n 1 ${runname}_Nextclade.results2.csv && tail -n +2 ${runname}_Nextclade.results2.csv | sort) > ${runname}_Nextclade.results2_sorted.csv
#(head -n 2 ${runname}_Nextclade.results2.csv && tail -n +3 ${runname}_Nextclade.results2.csv | sort) > ${runname}_Nextclade.results2_sorted.csv

paste -d, ${runname}_Nextclade.results2_sorted.csv pangolin_out_sorted.csv > NextcladeAndPangolin.out.csv
paste -d, NextcladeAndPangolin.out.csv nextclade.out2_sorted.csv > NextcladeAndPangolin.out2.csv

sed 's/,/\t/g' NextcladeAndPangolin.out2.csv | sed 's/ORF10/ORF10\t/g' > ${runname}_NextcladeAndPangolin.csv

rm *results*
rm *out*

echo "Pangolin og Nextclade ferdig"

####### Pangolin og Nextclade  ##### SLUTT #####

####### Ekstrahere Spike fra fasta-sekvens #######	START ######
#input spike ref og fasta-fil

cd "${basedir}/${runname}_summaries/fasta"

Rscript /home/docker/Scripts/CSAK_Spike_Extractor_docker.R "${runname}.fa"

####### Ekstrahere Spike fra fasta-sekvens #######	SLUTT ######

##QC-plot/Noise_extractor/csv_merger
cd "${basedir}/${runname}_summaries/"

Rscript /home/docker/Scripts/CSAK_csv_merger_Illumina_docker.R
Rscript /home/docker/Scripts/CSAK_NoiseExtractor_docker.R c8

mkdir Frameshift
cp ${basedir}/${runname}_summaries/fasta/${runname}.fa Frameshift
cd Frameshift
Rscript /home/docker/Scripts/CSAK_Frameshift_Finder_docker.R c8
mv *.xlsx ${basedir}/${runname}_summaries/
cd ..
rm -rf Frameshift

cd "${basedir}"
Rscript /home/docker/Scripts/CSAK_QCPlotter_docker.R

Rscript /home/docker/Scripts/CoronaTree.R
mv *aligned.fasta ${basedir}/${runname}_summaries/

mv /home/docker/Fastq/Tree.pdf ${basedir}/${runname}_summaries/${runname}_tree.pdf
Rscript /home/docker/Scripts/CoverageCalculator.R


rm /home/docker/Fastq/*.fasta
rm /home/docker/Fastq/*.csv
rm /home/docker/Fastq/primers.bed
rm -r /home/docker/Fastq/temp/

##
rm $CoronaRef
rm $PrimerBed
rm ${CoronaRef}.fai
rm /home/docker/Fastq/Rplots.pdf
