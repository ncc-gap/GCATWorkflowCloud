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
    -L=${INTERVAL_AUTOSOME} \
    -I=${INPUT_CRAM} \
    -O=${OUTPUT_DIR}/${SAMPLE}.autosome.g.vcf \
    -R=${REFERENCE_DIR}/${REFERENCE_FASTA} \
    --native-pair-hmm-threads=$(nproc) \
    --sample-ploidy 2 \
    --tmp-dir=${OUTPUT_DIR}_tmp

bgzip -@ $(nproc) ${OUTPUT_DIR}/${SAMPLE}.autosome.g.vcf
tabix -p vcf ${OUTPUT_DIR}/${SAMPLE}.autosome.g.vcf.gz
rm -f ${OUTPUT_DIR}/${SAMPLE}.autosome.g.vcf.idx

/usr/bin/java \
    -XX:-UseContainerSupport \
    -Xmx30g \
    -jar ${GATK_JAR} HaplotypeCaller \
    -ERC GVCF \
    -L=${INTERVAL_PAR} \
    -I=${INPUT_CRAM} \
    -O=${OUTPUT_DIR}/${SAMPLE}.PAR.g.vcf \
    -R=${REFERENCE_DIR}/${REFERENCE_FASTA} \
    --native-pair-hmm-threads=$(nproc) \
    --sample-ploidy 2 \
    --tmp-dir=${OUTPUT_DIR}_tmp

bgzip -@ $(nproc) ${OUTPUT_DIR}/${SAMPLE}.PAR.g.vcf
tabix -p vcf ${OUTPUT_DIR}/${SAMPLE}.PAR.g.vcf.gz
rm -f ${OUTPUT_DIR}/${SAMPLE}.PAR.g.vcf.idx

/usr/bin/java \
    -XX:-UseContainerSupport \
    -Xmx30g \
    -jar ${GATK_JAR} HaplotypeCaller \
    -ERC GVCF \
    -L=${INTERVAL_CHRX} \
    -I=${INPUT_CRAM} \
    -O=${OUTPUT_DIR}/${SAMPLE}.chrX.female.g.vcf \
    -R=${REFERENCE_DIR}/${REFERENCE_FASTA} \
    --native-pair-hmm-threads=$(nproc) \
    --sample-ploidy 2 \
    --tmp-dir=${OUTPUT_DIR}_tmp

bgzip -@ $(nproc) ${OUTPUT_DIR}/${SAMPLE}.chrX.female.g.vcf
tabix -p vcf ${OUTPUT_DIR}/${SAMPLE}.chrX.female.g.vcf.gz
rm -f ${OUTPUT_DIR}/${SAMPLE}.chrX.female.g.vcf.idx

/usr/bin/java \
    -XX:-UseContainerSupport \
    -Xmx30g \
    -jar ${GATK_JAR} HaplotypeCaller \
    -ERC GVCF \
    -L=${INTERVAL_CHRX} \
    -I=${INPUT_CRAM} \
    -O=${OUTPUT_DIR}/${SAMPLE}.chrX.male.g.vcf \
    -R=${REFERENCE_DIR}/${REFERENCE_FASTA} \
    --native-pair-hmm-threads=$(nproc) \
    --sample-ploidy 1 \
    --tmp-dir=${OUTPUT_DIR}_tmp

bgzip -@ $(nproc) ${OUTPUT_DIR}/${SAMPLE}.chrX.male.g.vcf
tabix -p vcf ${OUTPUT_DIR}/${SAMPLE}.chrX.male.g.vcf.gz
rm -f ${OUTPUT_DIR}/${SAMPLE}.chrX.male.g.vcf.idx

/usr/bin/java \
    -XX:-UseContainerSupport \
    -Xmx30g \
    -jar ${GATK_JAR} HaplotypeCaller \
    -ERC GVCF \
    -L=${INTERVAL_CHRY} \
    -I=${INPUT_CRAM} \
    -O=${OUTPUT_DIR}/${SAMPLE}.chrY.male.g.vcf \
    -R=${REFERENCE_DIR}/${REFERENCE_FASTA} \
    --native-pair-hmm-threads=$(nproc) \
    --sample-ploidy 1 \
    --tmp-dir=${OUTPUT_DIR}_tmp

bgzip -@ $(nproc) ${OUTPUT_DIR}/${SAMPLE}.chrY.male.g.vcf
tabix -p vcf ${OUTPUT_DIR}/${SAMPLE}.chrY.male.g.vcf.gz
rm -f ${OUTPUT_DIR}/${SAMPLE}.chrY.male.g.vcf.idx

