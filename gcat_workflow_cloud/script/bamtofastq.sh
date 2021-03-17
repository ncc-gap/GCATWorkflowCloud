#!/bin/bash

set -o errexit
set -o nounset
set -o pipefail
set -x

F1_NAME=${OUTPUT_DIR}"/1_1.fastq"
F2_NAME=${OUTPUT_DIR}"/1_2.fastq"
O1_NAME=${OUTPUT_DIR}'/unmatched_first_output.txt'
O2_NAME=${OUTPUT_DIR}'/unmatched_second_output.txt'
T=${OUTPUT_DIR}'/temp.txt'
S=${OUTPUT_DIR}'/single_end_output.txt'

mkdir -p ${OUTPUT_DIR}

export LD_LIBRARY_PATH=/usr/local/lib

for NUM in `seq 0 ${SAMPLE_MAX_INDEX}`
do
    INPUT_BAM="INPUT_BAM_${NUM}"

    touch ${F1_NAME} ${F2_NAME} ${T} ${S} ${O1_NAME} ${O2_NAME}

    /usr/local/bin/bamtofastq ${PARAM} filename=${!INPUT_BAM} F=${F1_NAME}.tmp F2=${F2_NAME}.tmp T=${T}.tmp S=${S}.tmp O=${O1_NAME}.tmp O2=${O2_NAME}.tmp
    if [ -s ${F1_NAME}.tmp ]; then
        cat ${F1_NAME}.tmp >> ${F1_NAME}
        cat ${F2_NAME}.tmp >> ${F2_NAME}
        if [ -s ${S}.tmp ]; then
            cat ${S}.tmp >> ${S}
        fi
    elif [ -s ${S}.tmp ]; then
        cat ${S}.tmp >> ${F1_NAME}
    fi
    if [ -s ${T}.tmp ]; then        
        cat ${T}.tmp >> ${T}
    fi
    if [ -s ${O1_NAME}.tmp ]; then
        cat ${O1_NAME}.tmp >> ${O1_NAME}
    fi
    if [ -s ${O2_NAME}.tmp ]; then
        cat ${O2_NAME}.tmp >> ${O2_NAME}
    fi

    rm -f ${F1_NAME}.tmp ${F2_NAME}.tmp ${T}.tmp ${S}.tmp ${O1_NAME}.tmp ${O2_NAME}.tmp
done
