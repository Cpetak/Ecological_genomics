#!/bin/bash

cd ~/Ecological_genomics/myresults/

mkdir fastqc_trimmed

for file in /data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/XCV*.cl.pd.fq

do

# fastqc is a program we are calling that is already available, followed by input and output location
fastqc ${file} -o fastqc_trimmed/ # because I'm going to be in myresults already

done



 
