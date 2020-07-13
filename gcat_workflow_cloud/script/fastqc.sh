#!/bin/bash

set -xv

set -o errexit
set -o nounset
set -o pipefail

mkdir -p ${OUTPUT_DIR}
/tools/FastQC/fastqc -t $(nproc) ${FASTQC_PARAMS} -o ${OUTPUT_DIR} ${INPUT_FASTQ_R1} ${INPUT_FASTQ_R2}
