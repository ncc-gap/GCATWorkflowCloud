#!/bin/bash

set -xv

set -o errexit
set -o nounset
set -o pipefail

if [ ! -d ${OUTPUT_DIR} ]; then
    mkdir -p ${OUTPUT_DIR}
fi

/usr/bin/java \
    -XX:-UseContainerSupport \
    -jar /tools/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar CollectWgsMetrics \
    -I=${INPUT_CRAM} \
    -O=${OUTPUT_DIR}/${SAMPLE_NAME}.CollectWgsMetrics.txt \
    -R=${REFERENCE_DIR}/${REFERENCE_FASTA} \
    --TMP_DIR=${OUTPUT_DIR}
