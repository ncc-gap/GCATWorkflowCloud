#!/bin/bash

set -xv

set -o errexit
set -o nounset
set -o pipefail

sample_name=${1}
input_fastq_1=${2}
input_fastq_2=${3}
reference_fasta=${4}
output_dir=${5}
tmp_dir=${6}
num_gpus=${7}

/opt/pkg/parabricks/v2.5.0/pbrun fq2bam \
    --ref ${reference_fasta} \
    --in-fq ${input_fastq_1} ${input_fastq_2} "@RG\tID:${sample_name}\tPL:na\tLB:na\tSM:${sample_name}\tPU:${sample_name}" \
    --out-bam ${output_dir}/${sample_name}.pbrun-fq2bam.bam \
    --bwa-options "-T 0" \
    --tmp-dir ${tmp_dir} \
    --num-gpus ${num_gpus}
