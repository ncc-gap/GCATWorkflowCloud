#!/bin/bash

set -o errexit
set -o nounset

mkdir -p ${OUTPUT_DIR}
export INPUT_BAM=${OUTPUT_DIR}/temp.bam

/tools/samtools-1.9/samtools view \
    -T ${REFERENCE_DIR}/${REFERENCE_FILE} \
    -h \
    -b \
    -@ $(nproc) \
    ${INPUT_CRAM} > ${INPUT_BAM}

/tools/samtools-1.9/samtools index \
    ${INPUT_BAM}

cd /opt/gridss

bash gridss.sh \
    -o ${OUTPUT_DIR}/${VCF}  \
    -a ${OUTPUT_DIR}/${ASSEMBLE} \
    -r ${REFERENCE_DIR}/${REFERENCE_FILE}  \
    -j ${GRIDSS_JAR} \
    -t $(nproc) \
    -w ${OUTPUT_DIR} \
    --picardoptions VALIDATION_STRINGENCY=LENIENT \
    ${INPUT_BAM}

rm -f ${INPUT_BAM}
rm -f ${INPUT_BAM}.bai
