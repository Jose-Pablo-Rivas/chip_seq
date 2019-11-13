## Autores: Martin Moreno-Perez y Jose Pablo Rivas-Fernandez
##
## Date: November 2019
##
## Contact: marmorper20@alum.us.es / josrivfer1@alum.us.es
##
## Bioinformatica y analisis genomico
##
## Analisis de datos de Chip-Seq

#! /bin/bash

# Reading input parameters:

INPUT_ID=$1
WD=$2
NI=$3
NS=$4
SCRIPT=$5

# Adress sample folder

cd $WD/samples/sample_input_${INPUT_ID}

## Sample quality control and read mapping to reference genome
fastqc input_${INPUT_ID}.fq.gz

bowtie -x $WD/genome/index -U input_${INPUT_ID}.fq.gz -S input_${INPUT_ID}.sam

## Generting sorted bam file

samtools view -@ 2 -S -b input_${INPUT_ID}.sam > input_${INPUT_ID}.bam
rm input_${INPUT_ID}.sam
samtools sort -o input_${INPUT_ID}.bam -o input_${INPUT_ID}_sorted.bam
samtools index input_${INPUT_ID}_sorted.bam


# Sincronizacion

echo "input_${INPUT_ID} DONE" >> $WD/logs/blackboard_chip_seq.txt

DONE_INPUT=$(cat $WD/logs/blackboard_chip_seq.txt | grep "DONE" | wc -l)

if [$DONE_INPUT -eq $NS]
then
   qsub -N callpeak -o $WD/logs/callpeak $SCRIPT/callpeak.sh $WD $NI
fi

