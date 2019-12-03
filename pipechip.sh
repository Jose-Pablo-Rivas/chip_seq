## Authors: Jose Pablo Rivas-Fernandez y Martin Moreno-Perez
##
## Date: November 2019
##
## Contact: josrivfer1@alum.us.es /  marmorper20@alum.us.es
##
## ChIP-Seq Data Analysis
##
## Pipeline script


#! /bin/bash

if [ $# -eq 0 ]
  then
   echo "This pipeline analyzes ChIP-seq Data"
   echo "Usage: pipechip <param_files>"
   echo ""
   echo "param_file: file with the parameters specifications. Please, check test/params.txt for an example" # OJO con donde esta params.txt
   echo ""
   echo "Have fun!"

   exit 0
fi

# Loading parameters

PARAMS=$1

WD=$(grep working_directory: $PARAMS | awk '{ print $2 }' )
NS=$(grep number_of_samples: $PARAMS | awk '{ print $2 }' )
QGNM=$(grep genome_question: $PARAMS | awk '{ print $2 }' )
GNM=$(grep genome: $PARAMS | awk '{ print $2 }' )
QANT=$(grep annotation_question: $PARAMS | awk '{ print $2 }' )
ANT=$(grep annotation: $PARAMS | awk '{ print $2 }' )
NC=$(grep number_of_chip: $PARAMS | awk '{ print $2 }' )
NI=$(grep number_of_input: $PARAMS | awk '{ print $2 }' )
SCRIPT=$(grep script: $PARAMS | awk '{ print $2 }' )
TEST=$(grep test: $PARAMS | awk '{ print $2 }' )



echo "Is this a test?" $TEST

echo "Working directory is in:" $WD
echo "Number of samples:" $NS

if [ $QGNM == "YES" ]
then
   echo "Genome copied from:" $GNM
else
   echo "Genome downloaded from" $GNM
fi


if [ $QANT == "YES" ]
then
   echo "Annotation copied from:" $ANT
else
   echo "Annotation downloaded from" $ANT
fi

echo "Number of chip samples:" $NC
echo "Number of input samples:" $NI
echo "Scripts are in:" $SCRIPT



#echo "El directorio para el espacio de trabajo es:" $WD
#echo "El numero de muestras es:" $NS
#echo "El genoma de referencia se ha obtenido de la URL:" $GNM
#echo "La anotacion se ha obtenido de la siguiente URL:" $ANT
#echo "El numero de muestras Chip es de:" $NC
#echo "El numero de muestras Input es de:" $NI
#echo "Los scripts están en:" $SCRIPT
#echo "¿Se trata de un test?" $TEST


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
      echo sample_chip_$(($J+1)) = ${SAMPLES[$I]}
      ((J++))
      ((I++))
   elif [ $K -lt $NI ]
   then
      echo sample_input_$(($K+1)) = ${SAMPLES[$I]}
      ((K++))
      ((I++))
   fi
done


# Generating working directory

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
      mkdir sample_chip_$J
      ((I++))
      ((J++))
   elif [ $K -le $NI ]
   then
      mkdir sample_input_$K
      ((I++))
      ((K++))
   fi
done


# Downloading/Copying samples

cd $WD/samples

if [ $TEST == "YES" ]
then
   I=0
   J=0
   K=0
   while [ $I -lt $NS ]
   do
      if [ $J -lt $NC ]
      then
         cp ${SAMPLES[$I]} sample_chip_$(($J + 1))/chip_$(($J + 1)).fastq
         echo chip_$(($J+1)) "copied"
         ((I++))
         ((J++))
      elif [ $K -lt $NI ]
      then
         cp ${SAMPLES[$I]} sample_input_$(($K + 1))/input_$(($K + 1)).fastq
         echo input_$(($K+1)) "copied"
         ((I++))
         ((K++))
      fi
   done
else
   I=0
   J=0
   K=0
   while [ $I -lt $NS ]
   do
      if [ $J -lt $NC ]
      then
         fastq-dump --split-files ${SAMPLES[$I]} -O ./sample_chip_$(($J+1))
         mv ./sample_chip_$(($J+1))/${SAMPLES[$I]}* ./sample_chip_$(($J+1))/chip_$(($J+1)).fastq
         echo chip_$(($J+1)) "downloaded"
         ((I++))
         ((J++))
      elif [ $K -lt $NI ]
      then
         fastq-dump --split-files ${SAMPLES[$I]} -O ./sample_input_$(($K+1))
         mv ./sample_input_$(($K+1))/${SAMPLES[$I]}* ./sample_input_$(($K+1))/input_$(($K+1)).fastq
         echo input_$(($K+1)) "downloaded"
         ((I++))
         ((K++))
      fi
   done
fi


echo "Samples labelled with 1 are the control ones"
echo "Samples labelled with 2 are the treatment  ones"


# Downloading reference genome

cd $WD/genome

if [ $QGNM == "YES" ]
then
   cp $GNM genome.fa
## We will only use the command "gunzip genome.fa.gz" if the file is packed.
   echo "Genome copied"
else
   wget -O genome.fa.gz $GNM
   gunzip genome.fa.gz
   echo "Genome downloaded"
fi


# Downloading annotation

cd $WD/annotation

if [ $QANT == "YES" ]
then
   cp $ANT annotation.gtf
## we wil only use the command "gunzip annotation.fa.gz" if the file is packed.
   echo "Annotation copied"
else
   wget -O annotation.gtf.gz $ANT
   gunzip annotation.fa.gz
   echo "Annotation downloaded"
fi


# Building reference genome index

cd $WD/genome
bowtie2-build genome.fa index


# Executing sample processing scripts

I=1
J=1
K=1
while [ $I -le $NS ]
do
   if [ $J -le $NC ]
   then
      qsub -N chip_$J -o $WD/logs/chip_$J $SCRIPT/chip_seq_chip_processing.sh $J $WD $NC $NS $SCRIPT
      ((I++))
      ((J++))
   elif [ $K -le $NI ]
   then
      qsub -N input_$K -o $WD/logs/input_$K $SCRIPT/chip_seq_input_processing.sh $K $WD $NI $NS $SCRIPT
      ((I++))
      ((K++))
   fi
done
