## Authors: Jose Pablo Rivas-Fernandez y Martin Moreno-Perez
##
## Date: November 2019
##
## Contact: josrivfer1@alum.us.es /  marmorper20@alum.us.es
##
## ChIP-Seq Data Analysis
##
## Chip processing script


#! /bin/bash

# Reading input parameters

CHIP_ID=$1
WD=$2
NC=$3
NS=$4
SCRIPT=$5


# Adress sample folder

cd $WD/samples/sample_chip_${CHIP_ID}


# Sample quality control and read mapping to reference genome

fastqc chip_${CHIP_ID}.fastq
bowtie2 -x $WD/genome/index -U chip_${CHIP_ID}.fastq -S chip_${CHIP_ID}.sam


# Generting sorted bam file

samtools view -@ 2 -S -b chip_${CHIP_ID}.sam > chip_${CHIP_ID}.bam
rm chip_${CHIP_ID}.sam
samtools sort chip_${CHIP_ID}.bam -o chip_${CHIP_ID}_sorted.bam
samtools index chip_${CHIP_ID}_sorted.bam


# Synchronization

echo "chip_${CHIP_ID} DONE" >> $WD/logs/blackboard_chip_seq.txt

DONE_CHIP=$(cat $WD/logs/blackboard_chip_seq.txt | grep "DONE" | wc -l)

if [ $DONE_CHIP -eq $NS ]
then
   qsub -N callpeak -o $WD/logs/callpeak $SCRIPT/callpeak.sh $WD $NC
fi

