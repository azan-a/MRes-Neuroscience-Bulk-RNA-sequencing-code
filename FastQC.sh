#!/bin/bash

# Specify which folder has the FASTQ files.
input_dir="./FASTQ Files"

#Specify which folder you want to store the report.
output_dir="./FASTQ Files report"

#Specify where the FastQC script is
fastqc_path="./fastqc_v0.12.1/FastQC/fastqc"

#Print out a message stating that we have began FastQC analysis
echo "Starting FastQC analysis..."

#This is a "for loop" which lets us go through all the FASTQ files in the folder
#we have specified the script to look into earlier (The input_dir folder)
for file in "$input_dir"/*.fastq.gz; do
  # Run FastQC
  "$fastqc_path" -o "$output_dir" "$file"
done

#This prints out a message telling us the script has finish processing all
#the FASTQ files.
echo "FastQC analysis completed :D"

