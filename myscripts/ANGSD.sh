#!/bin/bash

cd /users/c/p/cpetak/Ecological_genomics/myresults/ANGSD

mypop=XCV
REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

ANGSD -b XCV_bam.list \
-ref ${REF} -anc ${REF} \
-out ${mypop}_outFold \
-nThreads 1 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 20 \
-minQ 20 \
-setMinDepth 3 \
-minInd 2 \
-setMinDepthInd 1 \
-setMaxDepthInd 17 \
-skipTriallelic 1 \
-GL 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doSaf 1 \
-pest ${mypop}_outFold.sfs \
-doThetas 1 \
-fold 1


thetaStat do_stat ${mypop}_outFold.thetas.idx
