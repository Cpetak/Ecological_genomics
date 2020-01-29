#!/bin/bash

cd ~/Ecological_genomics/myresults/

mkdir fastqc

for file in /data/project_data/RS_ExomeSeq/fastq/edge_fastq/XCV*fastq.gz

do

# fastqc is a program we are calling that is already available, followed by input and output location
fastqc ${file} -o fastqc/ # because I'm going to be in myresults already

done



 
