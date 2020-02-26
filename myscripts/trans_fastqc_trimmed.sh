#!/bin/bash

cd ~/Ecological_genomics/myresults/

mkdir fastqc_trans_trimmed

for file in /data/project_data/RS_RNASeq/fastq/cleanreads/ASC_06_C*.cl.fq

do

# fastqc is a program we are calling that is already available, followed by input and output location
fastqc ${file} -o fastqc_trans_trimmed/ # because I'm going to be in myresults already

done

for file2 in /data/project_data/RS_RNASeq/fastq/cleanreads/ASC_06_D*.cl.fq

do

# fastqc is a program we are calling that is already available, followed by input and output location
fastqc ${file2} -o fastqc_trans_trimmed/ # because I'm going to be in myresults already

done

 
