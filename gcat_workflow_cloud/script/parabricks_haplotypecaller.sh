#!/bin/bash

set -xv

set -o errexit
set -o nounset
set -o pipefail

sample_name=${1}
input_bam=${2}
reference_fasta=${3}
output_dir=${4}
tmp_dir=${5}
num_gpus=${6}

/opt/pkg/parabricks/v2.5.0/pbrun haplotypecaller \
    --ref ${reference_fasta} \
    --in-bam ${input_bam} \
    --out-variants ${output_dir}/${sample_name}.pbrun-haplotypecaller.vcf \
    --tmp-dir ${tmp_dir} \
    --num-gpus ${num_gpus}
