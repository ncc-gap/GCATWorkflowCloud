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

if [ ${NPROC} = "proc" ]; then
    NPROC=$(nproc)
elif [ ${NPROC} = "proc_half" ]; then
    NPROC=$(echo $(nproc) | awk '{printf "%d", $1/2+0.5}')
fi

/usr/bin/java \
    -XX:-UseContainerSupport ${HAPLOTYPE_JAVA_OPTION} \
    -jar ${GATK_JAR} HaplotypeCaller \
    -ERC GVCF \
    -L=${INTERVAL} \
    -I=${INPUT_CRAM} \
    -O=${OUTPUT_DIR}/${SAMPLE}.${TAG}.g.vcf \
    -R=${REFERENCE_DIR}/${REFERENCE_FASTA} \
    --native-pair-hmm-threads=${NPROC} ${HAPLOTYPE_OPTION} \
    --sample-ploidy ${PLOIDY} \
    --tmp-dir=${OUTPUT_DIR}_tmp

bgzip -@ $(nproc) ${OUTPUT_DIR}/${SAMPLE}.${TAG}.g.vcf
tabix -p vcf ${OUTPUT_DIR}/${SAMPLE}.${TAG}.g.vcf.gz
rm -f ${OUTPUT_DIR}/${SAMPLE}.${TAG}.g.vcf.idx

