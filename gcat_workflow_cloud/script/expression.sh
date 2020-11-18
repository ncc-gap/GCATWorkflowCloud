#!/bin/bash

set -x
set -o errexit
set -o nounset

OUTPUT_PREF=${OUTPUT_DIR}/${SAMPLE}
mkdir -p ${OUTPUT_DIR}

samtools view ${INPUT_CRAM} -b -@ 6 -T ${REFERENCE} -o ${OUTPUT_PREF}.temp.bam
rm ${INPUT_CRAM}

featureCounts -T 4 -p -a ${GTF} -O -B -C -o ${OUTPUT_PREF}.txt ${OUTPUT_PREF}.temp.bam
rm ${OUTPUT_PREF}.temp.bam

python /tools/simple_exp/proc_fc.py ${OUTPUT_PREF}.txt ${OUTPUT_PREF}.txt.summary ${GTF} > ${OUTPUT_PREF}.txt.fpkm
