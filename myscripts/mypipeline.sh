#!/bin/bash

#Path to my repo:
myrepo="/users/c/p/cpetak/Ecological_genomics"

#My population:
mypop="XCV"

#Directory to our cleaned and paired reads:
input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${mypop}"

#Directory to store the outputs of our mapping:
output="/data/project_data/RS_ExomeSeq/mapping"

#Run mapping.sh
source ./mapping.sh

#Run the post-processing steps
source ./process_bam.sh
