ChIP-Sequencing (ChIP-Seq) is a method that combines chromatin immunoprecipitation with massive DNA sequencing techniques in 
order to identify the binding sites of DNA-associated proteins.

The script pipechip.sh contains a pipeline that analyzes ChIP-Seq data.
Quizá explicar por encima lo que hace el script:

  Create working directory with:
 - Genome
 - Annotation
 - Samples
 - Results
 - Logs 

The scripts chip_seq_chip_processing.sh and chip_seq_input_processing.sh process all samples in parallel.
The script callpeak.sh .........

File params.txt contains all the parameters required to analyze the data.

Have fun!
