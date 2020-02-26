#!/bin/bash

cd ~/Ecological_genomics/myresults/

mkdir fastqc_trans

for file in /data/project_data/RS_RNASeq/fastq/ASC_06_C*.fastq.gz

do

# fastqc is a program we are calling that is already available, followed by input and output location
fastqc ${file} -o fastqc/ # because I'm going to be in myresults already

done

for file2 in /data/project_data/RS_RNASeq/fastq/ASC_06_D*.fastq.gz

do

# fastqc is a program we are calling that is already available, followed by input and output location
fastqc ${file2} -o fastqc/ # because I'm going to be in myresults already

done

 
