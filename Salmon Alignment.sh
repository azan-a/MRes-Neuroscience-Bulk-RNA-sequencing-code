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

# Specify file paths (ensure paths with spaces are properly quoted)
input_dir="./trimmed_fastq_files"
output_dir="./salmon_output/"
index_dir="./rattus_norvegicus_index/"
SALMON_PATH="/nobackup/leedsomics_tools/salmon-latest_linux_x86_64/bin/salmon"

threads=12

# Loop over each .gz file in input directory
for file in "$input_dir"/*.gz; do    
    # Run Salmon Quantification
    $SALMON_PATH quant -i "$index_dir" -l A -r "$file" -p $threads -o "$output_dir" --validateMappings --seqBias
done

