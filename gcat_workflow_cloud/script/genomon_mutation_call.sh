#!/bin/bash

set -o errexit
set -o nounset

mkdir -p ${OUTPUT_DIR}
OUTPUT_PREF=${OUTPUT_DIR}/${SAMPLE1}
REFERENCE=${REFERENCE_DIR}/${REFERENCE_FILE}
SAMTOOLS=samtools

# INPUT_BAM1: target
INPUT_BAM1=${INPUT_BAM_DIR1}/${INPUT_BAM1}

# print_header
mut_header=""

if [ _${SAMPLE2} = "_None" ]; then 

    # Fisher's Exact Test
    if [ "_${FISHER_INTERVAL_LIST}" != "_" ]; then
        FISHER_SINGLE_OPTION="${FISHER_SINGLE_OPTION} -L ${FISHER_INTERVAL_LIST} "
    fi
    fisher single -o ${OUTPUT_PREF}.fisher_mutations.txt --ref_fa ${REFERENCE} -1 ${INPUT_BAM1} --samtools_path ${SAMTOOLS} ${FISHER_SINGLE_OPTION} --samtools_params "${FISHER_SAMTOOLS_OPTION}"

    # Local realignment using blat. The candidate mutations are varidated.
    mutfilter realignment --target_mutation_file ${OUTPUT_PREF}.fisher_mutations.txt -1 ${INPUT_BAM1} --output ${OUTPUT_PREF}.realignment_mutations.txt --ref_genome ${REFERENCE}  ${REALIGNMENT_OPTION}

    # Annotation if the candidate is on the simplerepeat. 
    mutfilter simplerepeat --target_mutation_file ${OUTPUT_PREF}.realignment_mutations.txt --output ${OUTPUT_PREF}.simplerepeat_mutations.txt --simple_repeat_db ${ANNOTATION_DB}/simpleRepeat.bed.gz

    # header
    mut_header="Chr,Start,End,Ref,Alt,depth,variantNum,bases,A_C_G_T,misRate,strandRatio,10%_posterior_quantile,posterior_mean,90%_posterior_quantile,readPairNum,variantPairNum,otherPairNum,10%_posterior_quantile(realignment),posterior_mean(realignment),90%_posterior_quantile(realignment),simple_repeat_pos,simple_repeat_seq"
    print_header=`echo $mut_header | tr "," "\t"`
    echo "$print_header" > ${OUTPUT_PREF}.genomon_mutation.result.txt
    cat ${OUTPUT_PREF}.simplerepeat_mutations.txt >> ${OUTPUT_PREF}.genomon_mutation.result.txt

    mutil filter -i ${OUTPUT_PREF}.genomon_mutation.result.txt -o ${OUTPUT_PREF}.genomon_mutation.result.filt.txt ${FILTER_SINGLE_OPTION}

else
    # INPUT_BAM2: pair
    INPUT_BAM2=${INPUT_BAM_DIR2}/${INPUT_BAM2}

    # Fisher's Exact Test
    if [ "_${FISHER_INTERVAL_LIST}" != "_" ]; then
        FISHER_PAIR_OPTION="${FISHER_PAIR_OPTION} -L ${FISHER_INTERVAL_LIST} "
    fi
    fisher comparison -o ${OUTPUT_PREF}.fisher_mutations.txt --ref_fa ${REFERENCE} -1 ${INPUT_BAM1} -2 ${INPUT_BAM2} --samtools_path ${SAMTOOLS} ${FISHER_PAIR_OPTION} --samtools_params "${FISHER_SAMTOOLS_OPTION}"
    
    # Local realignment using blat. The candidate mutations are varidated.
    mutfilter realignment --target_mutation_file ${OUTPUT_PREF}.fisher_mutations.txt -1 ${INPUT_BAM1} -2 ${INPUT_BAM2} --output ${OUTPUT_PREF}.realignment_mutations.txt --ref_genome ${REFERENCE} ${REALIGNMENT_OPTION}
    
    # Annotation if the candidate is near Indel. 
    mutfilter indel --target_mutation_file ${OUTPUT_PREF}.realignment_mutations.txt -2 ${INPUT_BAM2} --output ${OUTPUT_PREF}.indel_mutations.txt --ref_genome ${REFERENCE} --samtools_path ${SAMTOOLS} ${INDEL_OPTION} --samtools_params "${INDEL_SAMTOOLS_OPTION}"
    
    # Annotation if the candidate is near the breakpoint. 
    mutfilter breakpoint --target_mutation_file ${OUTPUT_PREF}.indel_mutations.txt -2 ${INPUT_BAM2} --output ${OUTPUT_PREF}.breakpoint_mutations.txt --ref_genome ${REFERENCE} ${BREAKPOINT_OPTION}
    
    # Annotation if the candidate is on the simplerepeat. 
    mutfilter simplerepeat --target_mutation_file ${OUTPUT_PREF}.breakpoint_mutations.txt --output ${OUTPUT_PREF}.simplerepeat_mutations.txt --simple_repeat_db ${ANNOTATION_DB}/simpleRepeat.bed.gz 

    # header
    mut_header="Chr,Start,End,Ref,Alt,depth_tumor,variantNum_tumor,depth_normal,variantNum_normal,bases_tumor,bases_normal,A_C_G_T_tumor,A_C_G_T_normal,misRate_tumor,strandRatio_tumor,misRate_normal,strandRatio_normal,P-value(fisher),readPairNum_tumor,variantPairNum_tumor,otherPairNum_tumor,readPairNum_normal,variantPairNum_normal,otherPairNum_normal,P-value(fisher_realignment),indel_mismatch_count,indel_mismatch_rate,bp_mismatch_count,distance_from_breakpoint,simple_repeat_pos,simple_repeat_seq"
    print_header=`echo $mut_header | tr "," "\t"`
    echo "$print_header" > ${OUTPUT_PREF}.genomon_mutation.result.txt
    cat ${OUTPUT_PREF}.simplerepeat_mutations.txt >> ${OUTPUT_PREF}.genomon_mutation.result.txt

    mutil filter -i ${OUTPUT_PREF}.genomon_mutation.result.txt -o ${OUTPUT_PREF}.genomon_mutation.result.filt.txt ${FILTER_PAIR_OPTION}
fi 

rm -f ${OUTPUT_PREF}.fisher_mutations.txt
rm -f ${OUTPUT_PREF}.realignment_mutations.txt
rm -f ${OUTPUT_PREF}.indel_mutations.txt
rm -f ${OUTPUT_PREF}.breakpoint_mutations.txt
rm -f ${OUTPUT_PREF}.simplerepeat_mutations.txt
