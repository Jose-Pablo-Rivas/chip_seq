ChIP-Sequencing (ChIP-Seq) is a method that combines chromatin immunoprecipitation with massive DNA sequencing techniques in 
order to identify the binding sites of DNA-associated proteins.

After downloading the material from GitHub, you must run the file - install.sh - .

This script is not executable. You can make it executable by typing the following command in the terminal and in the directory where 
the script is: " chmod +x install.sh ". After that, you will be able to run the script and follow the next steps. 

Install.sh will create a directory in which all the scripts and files will be saved. The directory will be created in your Desktop with the
name "pipechip_josemartin".
 
pipechip_josemartin constains the scripts and another directory called "test" where there will be a .zip with the prr5 samples and a R script. 
  
The script pipechip.sh contains a pipeline that analyzes ChIP-Seq data and submits to the queu the scripts chip_seq_chip_processing.sh
and chip_seq_input_processing.sh. These two scripts process all samples in parallel. 

All scripts explains in a little way what they are doing in each step with a short definition above each command. It could be good to have a look to the 
scripts before running them.

The file "params.txt" contains all the parameters requiered to analyze the data. The parameters by default are for the prr5 samples analysis.
If you want to make a different study, you must change the parameters to those that agree with your study. If you add more parameters, you must add more 
parameters to the pipeline.

In the file "params.txt" the words that are between < > need to be changed by the user. 

We also have a tutorial on YouTube in which you can learn how to install and run the pipeline. The link is under:

https://youtu.be/6VlC_B-KmVQ

We hope that the pipeline will be usefull for you. 

Enjoy it and have fun!
