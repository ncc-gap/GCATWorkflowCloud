#! /bin/bash
set -eux

OUTPUT_DIR=`dirname $OUTPUT_FILE`
mkdir -p ${OUTPUT_DIR}

samtools view ${INPUT_CRAM} -b -@ 6 -T ${REFERENCE} -o ${OUTPUT_DIR}/temp.bam
rm ${INPUT_CRAM}
samtools index ${OUTPUT_DIR}/temp.bam

intron_retention_utils simple_count ${OUTPUT_DIR}/temp.bam ${OUTPUT_DIR}/temp.txt --pass_bam_filt --genome_id hg38
rm ${OUTPUT_DIR}/temp.bam

tail -n +2 ${OUTPUT_DIR}/temp.txt | bgzip -c > ${OUTPUT_FILE}
tabix -p vcf ${OUTPUT_FILE}

rm ${OUTPUT_DIR}/temp.txt
