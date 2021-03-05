#!/bin/bash

set -o errexit
set -o nounset

mkdir -p ${OUTPUT_DIR}

if [ "${NORMAL_BAM_FILE}" != "" ]; then
    python /manta/bin/configManta.py \
        --normalBam ${NORMAL_BAM_DIR}/${NORMAL_BAM_FILE} \
        --tumorBam ${TUMOR_BAM_DIR}/${TUMOR_BAM_FILE} \
        --referenceFasta ${REFERENCE} \
        --runDir ${OUTPUT_DIR} ${CONFIG_MANTA_OPTION}
else
    python /manta/bin/configManta.py \
        --tumorBam ${TUMOR_BAM_DIR}/${TUMOR_BAM_FILE} \
        --referenceFasta ${REFERENCE} \
        --runDir ${OUTPUT_DIR} ${CONFIG_MANTA_OPTION}
fi

python ${OUTPUT_DIR}/runWorkflow.py -j $(nproc) ${WORKFLOW_OPTION}
