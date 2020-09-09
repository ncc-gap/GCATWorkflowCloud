#!/bin/bash

set -xv

set -o errexit
set -o nounset
set -o pipefail

if [ ! -d ${OUTPUT_DIR} ]; then
    mkdir -p ${OUTPUT_DIR}
fi
if [ ! -d ${OUTPUT_DIR}_tmp ]; then
    mkdir -p ${OUTPUT_DIR}_tmp
fi

/usr/bin/java \
    -XX:-UseContainerSupport \
    -Xmx30g \
    -jar ${GATK_JAR} HaplotypeCaller \
    -ERC GVCF \
    -L=${INTERVAL} \
    -I=${INPUT_CRAM} \
    -O=${OUTPUT_DIR}/${SAMPLE}.${TAG}.g.vcf \
    -R=${REFERENCE_DIR}/${REFERENCE_FASTA} \
    --native-pair-hmm-threads=$(nproc) \
    --sample-ploidy ${PLOIDY} \
    --tmp-dir=${OUTPUT_DIR}_tmp

bgzip -@ $(nproc) ${OUTPUT_DIR}/${SAMPLE}.${TAG}.g.vcf
tabix -p vcf ${OUTPUT_DIR}/${SAMPLE}.${TAG}.g.vcf.gz
rm -f ${OUTPUT_DIR}/${SAMPLE}.${TAG}.g.vcf.idx
