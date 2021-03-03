#!/bin/bash

set -o errexit
set -o nounset

mkdir -p ${OUTPUT_DIR}

python /manta/bin/configManta.py \
    --bam ${NORMAL_BAM_DIR}/${NORMAL_BAM_FILE} \
    --referenceFasta ${REFERENCE_DIR}/${REFERENCE_FILE} \
    --runDir ${OUTPUT_DIR} ${CONFIG_MANTA_OPTION}

python ${OUTPUT_DIR}/runWorkflow.py ${WORKFLOW_OPTION}

