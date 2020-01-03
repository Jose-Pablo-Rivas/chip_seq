ChIP-Sequencing (ChIP-Seq) is a method that combines chromatin immunoprecipitation with massive DNA sequencing techniques in 
order to identify the binding sites of DNA-associated proteins.

After downloading the material from GitHub, you must run the file install.sh. This file create a directory and save the scripts and files into it, to make the
way to work easier. The directory will be in your Desktop with the name "pipechip_josemartin".

In that directory, there will be the scripts and a directory called "test" which will have a .zip with the prr5 samples and a R script. 
  
The script pipechip.sh contains a pipeline that analyzes ChIP-Seq data and it submits to the queu the scripts chip_seq_chip_processing.sh
and chip_seq_input_processing.sh. These two scripts process all samples in parallel. 

All scripts explains in a little way what they are doing in each step with a short definition above each command. It could be good to have a look to the 
scripts before running them.

The file params.txt contains all the parameters requiered to analyze the data. The parameters that are by default, they are for the prr5 samples analyzes.
If you want to make a different study, you will change the parameters to those that agree with your study. If you add more parameters, you will add more 
parameters to the pipeline.

We hope that the pipeline will be usefull for you. 

Enjoy it and have fun!
