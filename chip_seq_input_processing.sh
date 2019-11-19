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
   qsub -N callpeak -o $WD/logs/callpeak $SCRIPT/callpeak.sh $WD $NI
fi
