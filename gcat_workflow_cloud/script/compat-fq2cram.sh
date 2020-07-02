#!/bin/bash

set -xv

set -o errexit
set -o nounset
set -o pipefail

work_dir=$(dirname ${OUTPUT_CRAM})

/tools/bwa-0.7.15/bwa mem \
    -t $(nproc) \
    -K 10000000 \
    -T 0 \
    -R "@RG\tID:${SAMPLE_NAME}\tPL:na\tLB:na\tSM:${SAMPLE_NAME}\tPU:${SAMPLE_NAME}" \
    ${REFERENCE_DIR}/${REFERENCE_FASTA} \
    ${INPUT_FASTQ_1} \
    ${INPUT_FASTQ_2} \
| /usr/bin/java \
    -XX:-UseContainerSupport \
    -Xmx30g \
    -jar /tools/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar SortSam \
    --MAX_RECORDS_IN_RAM=5000000 \
    -I=/dev/stdin \
    -O=${work_dir}/${SAMPLE_NAME}.bam \
    --SORT_ORDER=coordinate \
    --TMP_DIR=${work_dir}

/usr/bin/java \
    -XX:-UseContainerSupport \
    -Xmx30g \
    -jar /tools/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar MarkDuplicates \
    -I=${work_dir}/${SAMPLE_NAME}.bam \
    -O=${work_dir}/${SAMPLE_NAME}.markdup.bam \
    -M=${OUTPUT_MARKDUP_METRICS} \
    --TMP_DIR=${work_dir}

rm -f ${work_dir}/${SAMPLE_NAME}.bam

/tools/samtools-1.9/samtools view \
    -@ $(nproc) \
    -C \
    -T ${REFERENCE_DIR}/${REFERENCE_FASTA} \
    -o ${OUTPUT_CRAM} \
    ${work_dir}/${SAMPLE_NAME}.markdup.bam

/tools/samtools-1.9/samtools index \
    -@ $(nproc) \
    ${OUTPUT_CRAM}

rm -f ${work_dir}/${SAMPLE_NAME}.markdup.bam
