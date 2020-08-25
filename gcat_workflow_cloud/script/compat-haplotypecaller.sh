#!/bin/bash

set -xv

set -o errexit
set -o nounset
set -o pipefail

OUTPUT_DIR=$(dirname ${OUTPUT_VCF})
if [ ! -d ${OUTPUT_DIR} ]; then
    mkdir -p ${OUTPUT_DIR}
fi

/usr/bin/java \
    -XX:-UseContainerSupport \
    -Xmx30g \
    -jar ${GATK_JAR} HaplotypeCaller \
    -ERC GVCF \
    -L=${INTERVALS} \
    -I=${INPUT_CRAM} \
    -O=${OUTPUT_VCF} \
    -R=${REFERENCE_DIR}/${REFERENCE_FASTA} \
    --native-pair-hmm-threads=$(nproc) \
    --tmp-dir=${OUTPUT_DIR}
