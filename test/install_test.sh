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

if [ -d $HOME/opt/ ]
then
   mkdir $HOME/opt/pipechip_josemartin
   mv ./pipechip.sh $HOME/opt/pipechip_josemartin/
   mv ./chip_seq_chip_processing.sh $HOME/opt/pipechip_josemartin
   mv ./chip_seq_input_processing.sh $HOME/opt/pipechip_josemartin
   mv ./callpeak.sh $HOME/opt/pipechip_josemartin
   mkdir $HOME/opt/pipechip_josemartin/test/
   mv ./prr5_samples.zip $HOME/opt/pipechip_josemartin/test/
   cd $HOME/opt/pipechip_josemartin/test/
   unzip prr5_samples.zip
   rm prr5_samples.zip
else
   mkdir $HOME/opt/
   mkdir $HOME/opt/pipechip_josemartin
   mv ./pipechip.sh $HOME/opt/pipechip_josemartin/
   mv ./chip_seq_chip_processing.sh $HOME/opt/pipechip_josemartin
   mv ./chip_seq_input_processing.sh $HOME/opt/pipechip_josemartin
   mv ./callpeak.sh $HOME/opt/pipechip_josemartin
   mkdir $HOME/opt/pipechip_josemartin/test/
   mv ./prr5_samples.zip $HOME/opt/pipechip_josemartin/test/
   cd $HOME/opt/pipechip_josemartin/test/
   unzip prr5_samples.zip
   rm prr5_samples.zip
fi


# mv ./install_test.sh $HOME/opt/pipechip_josemartin/

