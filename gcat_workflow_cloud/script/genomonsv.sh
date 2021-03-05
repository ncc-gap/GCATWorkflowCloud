#!/bin/bash

set -o errexit
set -o nounset

mkdir -p ${TUMOR_OUTPUT_DIR}

TUMOR_BAM=${TUMOR_BAM_DIR}/${TUMOR_BAM}
NORMAL_BAM=${NORMAL_BAM_DIR}/${NORMAL_BAM}

# GenomonSV parse tumor
GenomonSV \
    parse \
    ${TUMOR_BAM} \
    ${TUMOR_OUTPUT_DIR}/${TUMOR_SAMPLE} \
    --reference ${REFERENCE} \
    ${GENOMONSV_PARSE_OPTION}


if [ ! _${NORMAL_SAMPLE} = "_None" ]
then
    mkdir -p ${NORMAL_OUTPUT_DIR}
    # GenomonSV parse normal
    GenomonSV \
        parse \
        ${NORMAL_BAM} \
        ${NORMAL_OUTPUT_DIR}/${NORMAL_SAMPLE} \
        --reference ${REFERENCE} \
        ${GENOMONSV_PARSE_OPTION}
fi

# GenomonSV filt
argument="${TUMOR_BAM} ${TUMOR_OUTPUT_DIR}/${TUMOR_SAMPLE} ${REFERENCE}"

if [ ! _${CONTROL_PANEL} = "_None" ]
then
    argument="${argument} --non_matched_control_junction ${MERGED_JUNCTION_DIR}/${CONTROL_PANEL}"
fi

if [ ! _${NORMAL_SAMPLE} = "_None" ]
then
    argument="${argument} --matched_control_bam ${NORMAL_BAM}"

    if [ ! _${CONTROL_PANEL} = "_None" ]
    then
        argument="${argument} --matched_control_label ${NORMAL_SAMPLE}"
    fi
fi

argument="${argument} ${GENOMONSV_FILT_OPTION}"

GenomonSV filt ${argument}

sv_utils filter ${TUMOR_OUTPUT_DIR}/${TUMOR_SAMPLE}.genomonSV.result.txt ${TUMOR_OUTPUT_DIR}/${TUMOR_SAMPLE}.genomonSV.result.filt.txt ${SV_UTILS_FILT_OPTION}
