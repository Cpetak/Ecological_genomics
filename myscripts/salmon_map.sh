#!/bin/bash

output=/data/project_data/RS_RNASeq/salmon/cleanedreads

cd /data/project_data/RS_RNASeq/fastq/cleanreads/

for file in ASC_06_C*.cl.fq

do

salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_HC27_index -l A -r ${file} --validateMappings -o ${output}/${file}

done

for file2 in ASC_06_D*.cl.fq

do

salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_HC27_index -l A -r ${file2} --validateMappings -o ${output}/${file2}

done
