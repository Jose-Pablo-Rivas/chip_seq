## Author: Jose Pablo Rivas-Fernandez y Martin Moreno-Perez
## Date: November 2019
## email: josrivfer1@alum.us.es

#! /bin/bash


if [ -d $HOME/opt/ ]
then
   mkdir $HOME/opt/pipechip_josemartin
   mv ../pipechip.sh $HOME/opt/pipechip_josemartin/
##Hay que añadir el resto de scripts
   mkdir $HOME/opt/pipechip_josemartin/test/
   mv ./prr5_samples.zip $HOME/opt/pipechip_josemartin/test/
   cd $HOME/opt/pipechip_josemartin/test/
   gunzip prr5_samples.zip
else
   mkdir $HOME/opt/
   mkdir $HOME/opt/pipechip_josemartin
   mv ../pipechip.sh $HOME/opt/pipechip_josemartin/
##Hay que añadir el resto de scripts
   mkdir $HOME/opt/pipechip_josemartin/test/
   mv ./prr5_samples.zip $HOME/opt/pipechip_josemartin/test/
   cd $HOME/opt/pipechip_josemartin/test/
   gunzip prr5_samples.zip
fi
