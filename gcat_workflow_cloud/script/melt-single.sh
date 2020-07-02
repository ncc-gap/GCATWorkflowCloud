#!/bin/bash

set -xv
set -o errexit
set -o nounset

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

ls /MELTv2.2.0/me_refs/Hg38/*.zip > mei_list.txt

time java \
    -XX:-UseContainerSupport \
    -jar /MELTv2.2.0/MELT.jar Single \
    -a \
    -h ${REFERENCE_DIR}/${REFERENCE_FILE} \
    -w ${OUTPUT_DIR} \
    -t mei_list.txt \
    -n /MELTv2.2.0/add_bed_files/Hg38/Hg38.genes.bed \
    -bamfile ${INPUT_BAM} \
    -bowtie /tools/bowtie2-2.4.1-linux-x86_64/bowtie2
