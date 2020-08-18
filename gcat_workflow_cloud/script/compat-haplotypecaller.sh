#!/bin/bash

set -xv

set -o errexit
set -o nounset
set -o pipefail

/usr/bin/java \
    -XX:-UseContainerSupport \
    -Xmx30g \
    -jar ${GATK_JAR} HaplotypeCaller \
    -I=${INPUT_CRAM} \
    -O=${OUTPUT_VCF} \
    -R=${REFERENCE_DIR}/${REFERENCE_FASTA} \
    --native-pair-hmm-threads=$(nproc) \
    --TMP_DIR=$(dirname ${OUTPUT_VCF})
