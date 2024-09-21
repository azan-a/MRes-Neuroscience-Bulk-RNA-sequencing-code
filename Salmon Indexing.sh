#!/bin/bash

#Run with current environment (-V) and in the current directory (-cwd)
#$ -V -cwd

#Request some time- min 15 mins - max 48 hours
#$ -l h_rt=12:00:00

#Request some memory per core
#$ -pe smp 2
#$ -l h_vmem=32G

#Get email at start and end of the job
#$ -m be

#Tell the script where it can find the following files
TRANSCRIPTOME_PLUS_GENOME="./transcriptome_plus_genome.fa.gz"
DECOYS="./decoys.txt"
OUTPUT_INDEX="./rattus_norvegicus_index/"
SALMON_PATH="/nobackup/leedsomics_tools/salmon-latest_linux_x86_64/bin/salmon"

#Create an index for salmon
$SALMON_PATH index -t $TRANSCRIPTOME_PLUS_GENOME -d $DECOYS -p 2 -i $OUTPUT_INDEX
