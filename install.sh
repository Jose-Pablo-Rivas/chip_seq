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

if [ -d $HOME/pipechip_josemartin/ ]
then
   mv ./pipechip.sh $HOME/pipechip_josemartin/
   mv ./chip_seq_chip_processing.sh $HOME/pipechip_josemartin
   mv ./chip_seq_input_processing.sh $HOME/pipechip_josemartin
   mv ./params.txt $HOME/pipechip_josemartin
   mv ./README.md $HOME/pipechip_josemartin
   mkdir $HOME/pipechip_josemartin/test/
   mv ./test/prr5_samples.zip $HOME/pipechip_josemartin/test/
   mv ./test/chip_seq_prr5.R $HOME/pipechip_josemartin/test/
   rm -r ./test
   cd $HOME/pipechip_josemartin/
   chmod +x pipechip.sh chip_seq_chip_processing.sh chip_seq_input_processing.sh
   cd $HOME/pipechip_josemartin/test/
   unzip prr5_samples.zip
   rm prr5_samples.zip
else
   mkdir $HOME/pipechip_josemartin
   mv ./pipechip.sh $HOME/pipechip_josemartin/
   mv ./chip_seq_chip_processing.sh $HOME/pipechip_josemartin
   mv ./chip_seq_input_processing.sh $HOME/pipechip_josemartin
   mv ./params.txt $HOME/pipechip_josemartin
   mv ./README.md $HOME/pipechip_josemartin
   mkdir $HOME/pipechip_josemartin/test/
   mv ./test/prr5_samples.zip $HOME/pipechip_josemartin/test/
   mv ./test/chip_seq_prr5.R $HOME/pipechip_josemartin/test/
   rm -r ./test
   cd $HOME/pipechip_josemartin/
   chmod +x pipechip.sh chip_seq_chip_processing.sh chip_seq_input_processing.sh
   cd $HOME/pipechip_josemartin/test/
   unzip prr5_samples.zip
   rm prr5_samples.zip
fi

cd $HOME/pipechip_josemartin

SWD=$(pwd)

echo "script:" $SWD >> params.txt

echo "The installation has finished"

rm ./install_test.sh
