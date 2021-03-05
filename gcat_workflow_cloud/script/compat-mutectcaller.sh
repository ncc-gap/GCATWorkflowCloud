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

if [ "${INPUT_NORMAL_CRAM}" != "" ]; then
    /usr/bin/java \
        -XX:-UseContainerSupport ${MUTECT_JAVA_OPTION} \
        -jar ${GATK_JAR} Mutect2 \
        -I=${INPUT_TUMOR_CRAM} -tumor ${TUMOR_SAMPLE} \
        -I=${INPUT_NORMAL_CRAM} -normal ${NORMAL_SAMPLE} \
        -O=${OUTPUT_DIR}/${TUMOR_SAMPLE}.mutectcaller.vcf \
        -R=${REFERENCE} \
        --native-pair-hmm-threads=${NPROC} ${MUTECT_OPTION} \
        --tmp-dir=${OUTPUT_DIR}_tmp
else
    /usr/bin/java \
        -XX:-UseContainerSupport ${MUTECT_JAVA_OPTION} \
        -jar ${GATK_JAR} Mutect2 \
        -I=${INPUT_TUMOR_CRAM} -tumor ${TUMOR_SAMPLE} \
        -O=${OUTPUT_DIR}/${TUMOR_SAMPLE}.mutectcaller.vcf \
        -R=${REFERENCE} ${MUTECT_OPTION} \
        --native-pair-hmm-threads=${NPROC} \
        --tmp-dir=${OUTPUT_DIR}_tmp
fi

bgzip -@ $(nproc) ${OUTPUT_DIR}/${TUMOR_SAMPLE}.mutectcaller.vcf
tabix -p vcf ${OUTPUT_DIR}/${TUMOR_SAMPLE}.mutectcaller.vcf.gz
rm -f ${OUTPUT_DIR}/${TUMOR_SAMPLE}.mutectcaller.vcf.idx

