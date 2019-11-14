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


## Reading parameters

WD=$1
NC=$2

## Llama a picos (cambiara dependiendo de los inputs y mocks que tengamos):

cd $WD/results

I=1

while [ $I -le $NC ]
do
   macs callpeak -t $WD/samples/sample_chip_$I/chip_${I}_sorted.bam -c $WD/samples/sample_input_$I/input_${I}_sorted.bam -n 'NAME' --outdir . -f BAM
   ((I++))
done


mv ../../chip_* $WD/results
mv ../../input_* $WD/results

echo "El analisis de los datos ha terminado" >> $WD/logs/finish.txt
