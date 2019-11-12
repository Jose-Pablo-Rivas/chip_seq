## Autores: Martin Moreno-Perez y Jose Pablo Rivas-Fernandez
##
## Date: November 2019
##
## Contact: marmorper20@alum.us.es / josrivfer1@alum.us.es
##
## Bioinformatica y analisis genomicos
##
## Analisis de datos de Chip-Seq

#! /bin/bash

if [ $# -eq 0 ]
  then
   echo "This pipeline analysis RNA-seq data"
   echo "Usage: piperna <param_files>"
   echo ""
   echo "param_file: file with the parameters specifications. Please, check test/params.txt for an example"
   echo ""
   echo "enjoy"

   exit 0
fi

##Loading parameters

PARAMS=$1

WD=$(grep working_directory: $PARAMS | awk '{ print $2 }' )
NS=$(grep number_of_samples: $PARAMS | awk '{ print $2 }' )
GNM=$(grep genome: $PARAMS | awk '{ print $2 }' )
ANT=$(grep annotation: $PARAMS | awk '{ print $2 }' )
NC=$(grep number_of_chip: $PARAMS | awk '{ print $2 }' )
NI=$(grep number_of_input: $PARAMS | awk '{ print $2 }' )
echo "El directorio para el espacio de trabajo es:" $WD
echo "El numero de muestras es:" $NS
echo "El genoma de referencia se ha obtenido de la URL:" $GNM
echo "La anotacion se ha obtenido de la siguiente URL:" $ANT
echo "El numero de muestras Chip es de:" $NC
echo "El numero de muestras Input es de:" $NI

SAMPLES=()

I=0
J=0
K=0

while [ $I -lt $NS ]
do
   if [ $J -lt $NC ]
   then
      SAMPLES[$I]=$(grep sample_chip_$(($J + 1)): $PARAMS | awk '{ print $2 }')
      ((I++))
      ((J++))
   elif [ $K -lt $NI ]
   then
      SAMPLES[$I]=$(grep sample_input_$(($K + 1)): $PARAMS | awk '{ print $2 }')
      ((I++))
      ((K++))
   fi
done

I=0
J=0
K=0

while [ $I -lt $NS ]
do
   if [ $J -lt $NC ]
   then
      echo sample_chip$(($J+1)) = ${SAMPLES[$I]}
      ((J++))
      ((I++))
   elif [ $K -lt $NI ]
   then
      echo sample_input$(($K+1)) = ${SAMPLES[$I]}
      ((K++))
      ((I++))
   fi
done


## Generate working directory
mkdir $WD
cd $WD
mkdir genome annotation samples results logs
cd samples

I=1
J=1
K=1
while [ $I -le $NS ]
do
   if [ $J -le $NC ]
   then
      mkdir sample_chip$J
      ((I++))
      ((J++))
   elif [ $K -le $NI ]
   then
       mkdir sample_input$K
      ((I++))
      ((K++))
   fi
done



##Copy samples

cd $WD/samples

I=0
J=0
K=0
while [ $I -lt $NS ]
do
   if [ $J -lt $NC ]
   then
      cp ${SAMPLES[$I]} sample_chip$(($J + 1))/Chip_$(($J + 1)).fq.gz
      ((I++))
      ((J++))
   elif [ $K -lt $NI ]
   then
      cp ${SAMPLES[$I]} sample_input$(($K + 1))/Input_$(($K + 1)).fq.gz
      ((I++))
      ((K++))
   fi
done

echo 'Las muestras 1 se corresponden con el control y las 2 con el tratamiento'


## Downloading reference genome

cd $WD/genome
cp $GNM genome.fa.gz
gunzip genome.fa.gz

## Downloading annotation

cd $WD/annotation
cp $ANT annotation.gtf.gz
gunzip annotation.fa.gz

## Building reference genome index

cd $WD/genome
bowtie2-build genome.fa index



