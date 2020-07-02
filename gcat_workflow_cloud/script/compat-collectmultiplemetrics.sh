#!/bin/bash

set -xv

set -o errexit
set -o nounset
set -o pipefail

if [ ! -d ${OUTPUT_DIR} ]; then
    mkdir -p ${OUTPUT_DIR}
fi

/usr/bin/java \
    -XX:-UseContainerSupport \
    -jar /tools/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar CollectMultipleMetrics \
    -I=${INPUT_CRAM} \
    -O=${OUTPUT_DIR}/${SAMPLE_NAME}.CollectMultipleMetrics \
    -R=${REFERENCE_DIR}/${REFERENCE_FASTA} \
    --TMP_DIR=${OUTPUT_DIR} \
    --PROGRAM CollectAlignmentSummaryMetrics \
    --PROGRAM CollectInsertSizeMetrics \
    --PROGRAM QualityScoreDistribution \
    --PROGRAM MeanQualityByCycle \
    --PROGRAM CollectBaseDistributionByCycle \
    --PROGRAM CollectGcBiasMetrics \
    --PROGRAM CollectSequencingArtifactMetrics \
    --PROGRAM CollectQualityYieldMetrics
