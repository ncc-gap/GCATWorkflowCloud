#!/bin/bash

set -xv

set -o errexit
set -o nounset
set -o pipefail

sample_name=${1}
input_tumor_bam=${2}
input_normal_bam=${3}
reference_fasta=${4}
output_dir=${5}
tmp_dir=${6}
pbrun=${7}

${pbrun} mutectcaller \
  --ref ${reference_fasta} ${mutect_option} \
  --in-tumor-bam ${input_tumor_bam} ${input_normal_bam} \
  --tumor-name ${sample_name} \
  --out-vcf ${output_dir}/${sample_name}.pbrun-mutectcaller.vcf \
  --tmp-dir ${tmp_dir}
