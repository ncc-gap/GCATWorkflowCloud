#!/bin/bash

set -xv

set -o errexit
set -o nounset
set -o pipefail

sample_name=${1}
input_bam=${2}
reference_fasta=${3}
output_dir=${4}
num_threads=${5}
bam_decompressor_threads=${6}

/opt/pkg/parabricks/v2.5.0/pbrun collectmultiplemetrics \
    --ref ${reference_fasta} \
    --bam ${input_bam} \
    --out-all-metrics ${output_dir}/${sample_name}.collectmultiplemetrics \
    --processor-threads ${num_threads} \
    --bam-decompressor-threads ${bam_decompressor_threads}
