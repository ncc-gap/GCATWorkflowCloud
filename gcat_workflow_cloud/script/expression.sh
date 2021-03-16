#!/bin/bash

set -x
set -o errexit
set -o nounset

OUTPUT_PREF=${OUTPUT_DIR}/${SAMPLE}
mkdir -p ${OUTPUT_DIR}

if [ ${SEQ_FORMAT} = "cram" ]; then
    BAM_FILE=${OUTPUT_PREF}.temp.bam
    samtools view ${INPUT_CRAM} -b -@ 6 -T ${REFERENCE} -o ${BAM_FILE}
    rm ${INPUT_CRAM}
else
    BAM_FILE=${INPUT_CRAM}
fi

featureCounts -T 4 -p -a ${GTF} -O -B -C -o ${OUTPUT_PREF}.txt ${BAM_FILE}
rm ${BAM_FILE}

python /tools/simple_exp/proc_fc.py ${OUTPUT_PREF}.txt ${OUTPUT_PREF}.txt.summary ${GTF} > ${OUTPUT_PREF}.txt.fpkm
