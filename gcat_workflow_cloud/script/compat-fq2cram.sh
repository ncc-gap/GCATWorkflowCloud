#!/bin/bash

set -xv

set -o errexit
set -o nounset
set -o pipefail

work_dir=$(dirname ${OUTPUT_CRAM})

array_rg=(${ARRAY_RG})

SORTED_BAMS=
REMOVE_BAMS=

for NUM in `seq 0 ${SAMPLE_MAX_INDEX}`
do
    INPUT_FASTQ_1="INPUT_FASTQ_1_${NUM}"
    INPUT_FASTQ_2="INPUT_FASTQ_2_${NUM}"
    RG=${array_rg[$NUM]}

    /tools/bwa-0.7.15/bwa mem \
        -t $(nproc) \
        -K 10000000 \
        -T 0 \
        -Y \
        -R "${RG}" \
        ${REFERENCE_DIR}/${REFERENCE_FASTA} \
        ${!INPUT_FASTQ_1} \
        ${!INPUT_FASTQ_2} \
    | /usr/bin/java \
        -XX:-UseContainerSupport \
        -Xmx30g \
        -jar ${GATK_JAR} SortSam \
        --MAX_RECORDS_IN_RAM=5000000 \
        -I=/dev/stdin \
        -O=${work_dir}/${SAMPLE_NAME}_${NUM}.bam \
        --SORT_ORDER=coordinate \
        --TMP_DIR=${work_dir}

    SORTED_BAMS=$SORTED_BAMS" -I=${work_dir}/${SAMPLE_NAME}_${NUM}.bam"
    REMOVE_BAMS=$REMOVE_BAMS" ${work_dir}/${SAMPLE_NAME}_${NUM}.bam"
    rm -f ${!INPUT_FASTQ_1} ${!INPUT_FASTQ_2}
done

/usr/bin/java \
    -XX:-UseContainerSupport \
    -Xmx30g \
    -jar ${GATK_JAR} MarkDuplicates \
    ${SORTED_BAMS} \
    -O=${work_dir}/${SAMPLE_NAME}.markdup.bam \
    -M=${OUTPUT_MARKDUP_METRICS} \
    --TMP_DIR=${work_dir}

rm -f ${REMOVE_BAMS}

if [ ${SEQ_FORMAT} = "cram" ]; then
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

else
    # index
    /usr/local/bin/samtools index \
        -@ $(nproc) \
        ${work_dir}/${SAMPLE_NAME}.markdup.bam
fi
