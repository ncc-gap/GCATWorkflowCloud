#!/bin/bash

set -xv

set -o errexit
set -o nounset
set -o pipefail

mkdir -p ${OUTPUT_DIR}

for i in `seq 0 ${SAMPLE_MAX_INDEX}`
do
    INPUT_FASTQ_1="INPUT_FASTQ_1_${i}"
    INPUT_FASTQ_2="INPUT_FASTQ_2_${i}"

    /tools/FastQC/fastqc -t $(nproc) ${FASTQC_PARAMS} -o ${OUTPUT_DIR} ${!INPUT_FASTQ_1} ${!INPUT_FASTQ_2}
done
