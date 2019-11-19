## Authors: Jose Pablo Rivas-Fernandez y Martin Moreno-Perez
##
## Date: November 2019
##
## Contact: josrivfer1@alum.us.es / marmorper20@alum.us.es
##
## ChIP-Seq Data Analysis
##
## Callpeak script


#! /bin/bash


# Reading input parameters

WD=$1
NC=$2


# Callpeak function

cd $WD/results

I=1

while [ $I -le $NC ]
do
   macs2 callpeak -t $WD/samples/sample_chip_$I/chip_${I}_sorted.bam -c $WD/samples/sample_input_$I/input_${I}_sorted.bam -n 'Picos_$I' --outdir . -f BAM
   ((I++))
done


# Cleaning XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

mv ../../chip_* $WD/results
mv ../../input_* $WD/results
mv ../../callpeak* $WD/results


# Finish message

echo "Data Analyis is finished!" >> $WD/logs/finish.txt
