## Authors: Jose Pablo Rivas-Fernandez y Martin Moreno-Perez
##
## Date: November 2019
##
## Contact: josrivfer1@alum.us.es /  marmorper20@alum.us.es
##
## ChIP-Seq Data Analysis
##
## Input procesing script


#! /bin/bash

# Reading input parameters

INPUT_ID=$1
WD=$2
NI=$3
NS=$4
SCRIPT=$5
LPROM=$6


# Adress sample folder

cd $WD/samples/sample_input_${INPUT_ID}


# Sample quality control and read mapping to reference genome

fastqc input_${INPUT_ID}.fastq

bowtie2 -x $WD/genome/index -U input_${INPUT_ID}.fastq -S input_${INPUT_ID}.sam


# Generting sorted bam file

samtools view -@ 2 -S -b input_${INPUT_ID}.sam > input_${INPUT_ID}.bam
rm input_${INPUT_ID}.sam
samtools sort input_${INPUT_ID}.bam -o input_${INPUT_ID}_sorted.bam
samtools index input_${INPUT_ID}_sorted.bam


# Synchronization

echo "input_${INPUT_ID} DONE" >> $WD/logs/blackboard_chip_seq.txt

DONE_INPUT=$(cat $WD/logs/blackboard_chip_seq.txt | grep "DONE" | wc -l)


if [ $DONE_INPUT -eq $NS ]
then
   # Callpeak function
   cd $WD/results
   I=1
   while [ $I -le $NC ]
   do
      macs2 callpeak -t $WD/samples/sample_chip_$I/chip_${I}_sorted.bam -c $WD/samples/sample_input_$I/input_${I}_sorted.bam -n Picos_$I --outdir . -f BAM
      ((I++))
   done
   # Cleaning 
   mv ../../chip_* $WD/results
   mv ../../input_* $WD/results
   # Finish message
   echo "Data Analyis is finished!" >> $WD/logs/finish.txt

   ## Analisys of the R Script

   mkdir $WD/results/R_analysis

   Rscript $SCRIPT/test/chip_seq_prr5.R $WD/results/Picos_1_peaks.narrowPeak $WD/results/R_analysis $LPROM

   ## Motifs analises by Hommer

   mkdir $WD/results/Homer
   findMotifsGenome.pl $WD/results/Picos_1_summits.bed tair10 $WD/results/Homer -size 60
fi


