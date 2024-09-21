#!/bin/bash

#Specify which folder has the FASTQ file.
input_dir="./FASTQ_files"

#Specify which folder you want the script to store the processed FASTQ file
output_dir="./Trimmed_FASTQ_files"

#Specify what adapter sequence you want to trim.
adapter_file="./Trimmomatic-0.39/adapters/NexteraPE-PE.fa"

#Specify where the trimmomatic executable file is
trimmomatic_jar="./Trimmomatic-0.39/trimmomatic-0.39.jar"

#Specify how many threads you want the software to use. Usually, you would set
#this value to the number of cores you have in your CPU. 
threads=12

#This is a "for loop" which lets us go through all the FASTQ files in the folder
#we have specified the script to look into earlier (The input_dir folder)
for file in "$input_dir"/*.gz; do

  #This extract the basename of the file, removing the .gz extension
  base_name=$(basename "$file" .gz)
  #Names the output file to basename + _trimmed.fastq.gz
  output_file="$output_dir/${base_name}_trimmed.fastq.gz"

  #Executes the trimmomatic command
  java -jar "$trimmomatic_jar" SE -threads $threads -phred33 "$file" \
  "$output_file" ILLUMINACLIP:"$adapter_file":2:30:8

  #Prints a message after each file has been processed
  echo "Processed $file"
done

#This prints out a message telling us the script has finish processing all
#the FASTQ files.
echo "All files processed :D"

