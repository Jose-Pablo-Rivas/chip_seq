## Authors: Jose Pablo Rivas-Fernandez y Martin Moreno-Perez
##
## Date: November 2019
##
## Contact: josrivfer1@alum.us.es /  marmorper20@alum.us.es
##
## ChIP-Seq Data Analysis
##
## Install script


#! /bin/bash

if [ -d $HOME/prueba/ ]
then
   mkdir $HOME/prueba/pipechip_josemartin
   mv ./pipechip.sh $HOME/prueba/pipechip_josemartin/
   mv ./chip_seq_chip_processing.sh $HOME/prueba/pipechip_josemartin
   mv ./chip_seq_input_processing.sh $HOME/prueba/pipechip_josemartin
   mv ./params.txt $HOME/prueba/pipechip_josemartin
   mv ./README.md $HOME/prueba/pipechip_josemartin
   mkdir $HOME/prueba/pipechip_josemartin/test/
   mv ./test/prr5_samples.zip $HOME/prueba/pipechip_josemartin/test/
   rm -r ./test
   cd $HOME/prueba/pipechip_josemartin/test/
   unzip prr5_samples.zip
   rm prr5_samples.zip
else
   mkdir $HOME/prueba/
   mkdir $HOME/prueba/pipechip_josemartin
   mv ./pipechip.sh $HOME/prueba/pipechip_josemartin/
   mv ./chip_seq_chip_processing.sh $HOME/prueba/pipechip_josemartin
   mv ./chip_seq_input_processing.sh $HOME/prueba/pipechip_josemartin
   mv ./params.txt $HOME/prueba/pipechip_josemartin
   mv ./README.md $HOME/prueba/pipechip_josemartin
   mkdir $HOME/prueba/pipechip_josemartin/test/
   mv ./test/prr5_samples.zip $HOME/prueba/pipechip_josemartin/test/
   rm -r ./test
   cd $HOME/prueba/pipechip_josemartin/test/
   unzip prr5_samples.zip
   rm prr5_samples.zip
fi

cd $HOME/prueba/pipechip_josemartin

SWD=$(pwd)

echo "script:" $SWD >> params.txt

echo "The installation has finished"

# rm ./install_test.sh
