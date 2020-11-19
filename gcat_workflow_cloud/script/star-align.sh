#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -o errexit
set -o pipefail
set -x

OUTPUT_PREF=${OUTPUT_DIR}/${SAMPLE}
mkdir -p ${OUTPUT_DIR}

# cat fastq
cat $(printenv | grep FASTQ1_ | egrep -v S3_ | cut -f 2 -d = | sort) > ${OUTPUT_DIR}/temp1.fastq
rm $(printenv | grep FASTQ1_ | egrep -v S3_ | cut -f 2 -d =)

# STAR
if [ ${FQ_TYPE} = "pair" ]; then
  cat $(printenv | grep FASTQ2_ | egrep -v S3_ | cut -f 2 -d = | sort) > ${OUTPUT_DIR}/temp2.fastq
  rm $(printenv | grep FASTQ2_ | egrep -v S3_ | cut -f 2 -d =)

  /usr/local/bin/STAR \
    --genomeDir ${STAR_REFERENCE} \
    --readFilesIn ${OUTPUT_DIR}/temp1.fastq ${OUTPUT_DIR}/temp2.fastq \
    --outFileNamePrefix ${OUTPUT_PREF}. \
    ${STAR_OPTION}

  rm ${OUTPUT_DIR}/temp2.fastq

else
  /usr/local/bin/STAR \
    --genomeDir ${STAR_REFERENCE} \
    --readFilesIn ${OUTPUT_DIR}/temp1.fastq \
    --outFileNamePrefix ${OUTPUT_PREF}. \
    ${STAR_OPTION}
fi

# remove fastq
rm ${OUTPUT_DIR}/temp1.fastq

# sort
/usr/local/bin/samtools sort \
  -T ${OUTPUT_PREF}.Aligned.sortedByCoord.out \
  -@ 6 -m 3G \
  ${OUTPUT_PREF}.Aligned.out.bam \
  -O bam > ${OUTPUT_PREF}.Aligned.sortedByCoord.out.bam

rm ${OUTPUT_PREF}.Aligned.out.bam

# cram
export REF_CACHE=/scratch/.cache/
mkdir -p ${REF_CACHE}

cd ${OUTPUT_DIR}
/usr/local/bin/samtools view \
    ${OUTPUT_PREF}.Aligned.sortedByCoord.out.bam \
    -@ 6 -C -T ${REFERENCE} \
    -o ${OUTPUT_PREF}.Aligned.sortedByCoord.out.cram

rm ${OUTPUT_PREF}.Aligned.sortedByCoord.out.bam

# index
/usr/local/bin/samtools index \
    ${OUTPUT_PREF}.Aligned.sortedByCoord.out.cram

