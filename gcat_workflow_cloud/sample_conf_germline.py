#! /usr/bin/env python
import gcat_workflow_cloud.sample_conf_abc as abc

class Sample_conf(abc.Sample_conf_abc):
    SECTION_FASTQ = "fastq"
    SECTION_BAM_IMPORT = "bam_import"
    SECTION_BAM_TOFASTQ = "bam_tofastq"
    SECTION_HTCALL = {
        "chr1": "haplotypecaller_parabricks_chr1",
        "chr2": "haplotypecaller_parabricks_chr2",
        "chr3": "haplotypecaller_parabricks_chr3",
        "chr4": "haplotypecaller_parabricks_chr4",
        "chr5": "haplotypecaller_parabricks_chr5",
        "chr6": "haplotypecaller_parabricks_chr6",
        "chr7": "haplotypecaller_parabricks_chr7",
        "chr8": "haplotypecaller_parabricks_chr8",
        "chr9": "haplotypecaller_parabricks_chr9",
        "chr10": "haplotypecaller_parabricks_chr10",
        "chr11": "haplotypecaller_parabricks_chr11",
        "chr12": "haplotypecaller_parabricks_chr12",
        "chr13": "haplotypecaller_parabricks_chr13",
        "chr14": "haplotypecaller_parabricks_chr14",
        "chr15": "haplotypecaller_parabricks_chr15",
        "chr16": "haplotypecaller_parabricks_chr16",
        "chr17": "haplotypecaller_parabricks_chr17",
        "chr18": "haplotypecaller_parabricks_chr18",
        "chr19": "haplotypecaller_parabricks_chr19",
        "chr20": "haplotypecaller_parabricks_chr20",
        "chr21": "haplotypecaller_parabricks_chr21",
        "chr22": "haplotypecaller_parabricks_chr22",
        "chrx_female": "haplotypecaller_parabricks_chrx_female",
        "chrx_male": "haplotypecaller_parabricks_chrx_male",
        "chry_male": "haplotypecaller_parabricks_chry_male",
        "par": "haplotypecaller_parabricks_par",
    }
    PAR = "haplotypecaller_parabricks_par"
    SECTION_WGS_METRICS = "collect_wgs_metrics"
    SECTION_MULTIPLE_METRICS = "collect_multiple_metrics"
    SECTION_GRIDSS = "gridss"
    SECTION_MANTA = "manta"
    SECTION_MELT = "melt"
    SECTION_READGROUP = "readgroup"
    SECTION_FASTQC = "fastqc"

    def __init__(self, sample_conf_file, exist_check = True):

        self.fastq = {}
        self.fastq_src = {}
        self.bam_tofastq = {}
        self.bam_tofastq_src = {}
        self.bam_import = {}
        self.bam_import_src = {}
        self.haplotype_call_chr1 = []
        self.haplotype_call_chr2 = []
        self.haplotype_call_chr3 = []
        self.haplotype_call_chr4 = []
        self.haplotype_call_chr5 = []
        self.haplotype_call_chr6 = []
        self.haplotype_call_chr7 = []
        self.haplotype_call_chr8 = []
        self.haplotype_call_chr9 = []
        self.haplotype_call_chr10 = []
        self.haplotype_call_chr11 = []
        self.haplotype_call_chr12 = []
        self.haplotype_call_chr13 = []
        self.haplotype_call_chr14 = []
        self.haplotype_call_chr15 = []
        self.haplotype_call_chr16 = []
        self.haplotype_call_chr17 = []
        self.haplotype_call_chr18 = []
        self.haplotype_call_chr19 = []
        self.haplotype_call_chr20 = []
        self.haplotype_call_chr21 = []
        self.haplotype_call_chr22 = []
        self.haplotype_call_chrx_female = []
        self.haplotype_call_chrx_male = []
        self.haplotype_call_chry_male = []
        self.haplotype_call_par = []
        self.wgs_metrics = []
        self.multiple_metrics = []
        self.gridss = []
        self.manta = []
        self.melt = []
        self.readgroup = {}
        self.readgroup_src = {}
        self.fastqc = []
        self.exist_check = exist_check
        
        self.parse_file(sample_conf_file)

    def parse_data(self, _data):
        
        input_sections = [self.SECTION_FASTQ, self.SECTION_BAM_IMPORT, self.SECTION_BAM_TOFASTQ]
        analysis_sections = [
            self.SECTION_WGS_METRICS, self.SECTION_MULTIPLE_METRICS, 
            self.SECTION_GRIDSS, self.SECTION_MANTA, self.SECTION_MELT,
            self.SECTION_FASTQC, 
        ]
        for key in self.SECTION_HTCALL:
           analysis_sections.append(self.SECTION_HTCALL[key])

        controlpanel_sections = []
        extend_sections = [self.SECTION_READGROUP]
        splited = self.split_section_data(_data, input_sections, analysis_sections, controlpanel_sections, extend_sections)
        
        bwa_samples = []
        if self.SECTION_FASTQ in splited:
            parsed_fastq = self.parse_data_fastq_pair(splited[self.SECTION_FASTQ])
            self.fastq.update(parsed_fastq["fastq"])
            self.fastq_src.update(parsed_fastq["fastq_src"])
            bwa_samples += self.fastq.keys()
            
        if self.SECTION_BAM_TOFASTQ in splited:
            parsed_bam_tofastq = self.parse_data_bam_tofastq(splited[self.SECTION_BAM_TOFASTQ])
            self.bam_tofastq.update(parsed_bam_tofastq["bam_tofastq"])
            self.bam_tofastq_src.update(parsed_bam_tofastq["bam_tofastq_src"])
            bwa_samples += self.bam_tofastq.keys()
            
        if self.SECTION_BAM_IMPORT in splited:
            parsed_bam_import = self.parse_data_bam_import(splited[self.SECTION_BAM_IMPORT])
            self.bam_import.update(parsed_bam_import["bam_import"])
            self.bam_import_src.update(parsed_bam_import["bam_import_src"])
            
        for key in self.SECTION_HTCALL:
            if self.SECTION_HTCALL[key] in splited:
                self.__dict__["haplotype_call_" + key] += self.parse_data_general(splited[self.SECTION_HTCALL[key]])

        if self.SECTION_WGS_METRICS in splited:
            self.wgs_metrics += self.parse_data_general(splited[self.SECTION_WGS_METRICS])

        if self.SECTION_MULTIPLE_METRICS in splited:
            self.multiple_metrics += self.parse_data_general(splited[self.SECTION_MULTIPLE_METRICS])

        if self.SECTION_GRIDSS in splited:
            self.gridss += self.parse_data_general(splited[self.SECTION_GRIDSS])

        if self.SECTION_MANTA in splited:
            self.manta += self.parse_data_general(splited[self.SECTION_MANTA])
        
        if self.SECTION_MELT in splited:
            self.melt += self.parse_data_general(splited[self.SECTION_MELT])

        if self.SECTION_FASTQC in splited:
            self.fastqc += self.parse_data_general(splited[self.SECTION_FASTQC])

        if len(bwa_samples) > 0:
            if not self.SECTION_READGROUP in splited:
                err_msg = "[%s] section is not set" % (self.SECTION_READGROUP)
                raise ValueError(err_msg)
            parsed_bam_tofastq = self.parse_data_readgroup(splited[self.SECTION_READGROUP], bwa_samples, self.SECTION_READGROUP)
            self.readgroup.update(parsed_bam_tofastq["metadata"])
            self.readgroup_src.update(parsed_bam_tofastq["metadata_src"])
