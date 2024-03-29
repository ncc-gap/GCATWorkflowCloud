#! /usr/bin/env python
import gcat_workflow_cloud.sample_conf_abc as abc

class Sample_conf(abc.Sample_conf_abc):
    SECTION_FASTQ = "fastq"
    SECTION_BAM_IMPORT = "bam_import"
    SECTION_BAM_TOFASTQ = "bam_tofastq"
    SECTION_MTCALL = "mutectcaller_parabricks"
    SECTION_WGS_METRICS = "collect_wgs_metrics"
    SECTION_MULTIPLE_METRICS = "collect_multiple_metrics"
    SECTION_GRIDSS = "gridss"
    SECTION_MANTA = "manta"
    SECTION_GENOMON_SV = "genomon_sv"
    SECTION_GENOMON_MUTATION_CALL = "genomon_mutation_call"
    SECTION_READGROUP = "readgroup"
    SECTION_FASTQC = "fastqc"

    def __init__(self, sample_conf_file, exist_check = True, use_bam = False):

        self.fastq = {}
        self.fastq_src = {}
        self.bam_tofastq = {}
        self.bam_tofastq_src = {}
        self.bam_import = {}
        self.bam_import_src = {}
        self.control_panel = {}
        self.mutect_call = []
        self.wgs_metrics = []
        self.multiple_metrics = []
        self.gridss = []
        self.manta = []
        self.genomon_sv = []
        self.genomon_mutation_call = []
        self.readgroup = {}
        self.readgroup_src = {}
        self.fastqc = []
        self.exist_check = exist_check
        self.use_bam = use_bam 
        self.parse_file(sample_conf_file)

    def parse_data(self, _data):
        
        input_sections = [self.SECTION_FASTQ, self.SECTION_BAM_IMPORT, self.SECTION_BAM_TOFASTQ]
        analysis_sections = [
            self.SECTION_MTCALL, 
            self.SECTION_WGS_METRICS, 
            self.SECTION_MULTIPLE_METRICS, 
            self.SECTION_GRIDSS, 
            self.SECTION_MANTA, 
            self.SECTION_GENOMON_SV,
            self.SECTION_GENOMON_MUTATION_CALL,
            self.SECTION_FASTQC, 
        ]
        controlpanel_sections = []
        extend_sections = [self.SECTION_READGROUP]
        splited = self.split_section_data(_data, input_sections, analysis_sections, controlpanel_sections, extend_sections)
        
        sample_ids = []
        bwa_samples = []
        if self.SECTION_FASTQ in splited:
            parsed_fastq = self.parse_data_fastq_pair(splited[self.SECTION_FASTQ])
            self.fastq.update(parsed_fastq["fastq"])
            self.fastq_src.update(parsed_fastq["fastq_src"])
            sample_ids += parsed_fastq["fastq"].keys()
            bwa_samples += self.fastq.keys()

        if self.SECTION_BAM_TOFASTQ in splited:
            parsed_bam_tofastq = self.parse_data_bam_tofastq(splited[self.SECTION_BAM_TOFASTQ])
            self.bam_tofastq.update(parsed_bam_tofastq["bam_tofastq"])
            self.bam_tofastq_src.update(parsed_bam_tofastq["bam_tofastq_src"])
            sample_ids += parsed_bam_tofastq["bam_tofastq"].keys()
            bwa_samples += self.bam_tofastq.keys()

        if self.SECTION_BAM_IMPORT in splited:
            parsed_bam_import = self.parse_data_bam_import(splited[self.SECTION_BAM_IMPORT])
            ext = "cram"
            if self.use_bam:
                ext = "bam"
            for sample in parsed_bam_import["bam_import"]:
                for path in parsed_bam_import["bam_import"][sample].split(";"):
                    if not path.rstrip().endswith(ext):
                        err_msg = "[%s]:%s, use %s file" % (self.SECTION_BAM_IMPORT, sample, ext)
                        raise ValueError(err_msg)
            self.bam_import.update(parsed_bam_import["bam_import"])
            self.bam_import_src.update(parsed_bam_import["bam_import_src"])
            sample_ids += parsed_bam_import["bam_import"].keys()
            
        if self.SECTION_MTCALL in splited:
            self.mutect_call += self.parse_data_tumor_normal(splited[self.SECTION_MTCALL], sample_ids, self.SECTION_MTCALL)
        
        if self.SECTION_WGS_METRICS in splited:
            self.wgs_metrics += self.parse_data_general(splited[self.SECTION_WGS_METRICS])

        if self.SECTION_MULTIPLE_METRICS in splited:
            self.multiple_metrics += self.parse_data_general(splited[self.SECTION_MULTIPLE_METRICS])

        if self.SECTION_GRIDSS in splited:
            self.gridss += self.parse_data_tumor_normal(splited[self.SECTION_GRIDSS], sample_ids, self.SECTION_GRIDSS)
            
        if self.SECTION_MANTA in splited:
            self.manta += self.parse_data_tumor_normal(splited[self.SECTION_MANTA], sample_ids, self.SECTION_MANTA)
            
        if self.SECTION_GENOMON_SV in splited:
            self.genomon_sv += self.parse_data_tumor_normal_controlpanel(
                splited[self.SECTION_GENOMON_SV],
                sample_ids,
                self.control_panel.keys(),
                self.SECTION_GENOMON_SV
            )

        if self.SECTION_GENOMON_MUTATION_CALL in splited:
            self.genomon_mutation_call += self.parse_data_tumor_normal_controlpanel(
                splited[self.SECTION_GENOMON_MUTATION_CALL],
                sample_ids,
                self.control_panel.keys(),
                self.SECTION_GENOMON_MUTATION_CALL
            )
        if self.SECTION_FASTQC in splited:
            self.fastqc += self.parse_data_general(splited[self.SECTION_FASTQC])
            for sample in self.fastqc:
                if sample in self.bam_import:
                    err_msg = "[%s] section, %s is not defined" % (self.SECTION_FASTQC, sample)
                    raise ValueError(err_msg)

        if len(bwa_samples) > 0:
            if not self.SECTION_READGROUP in splited:
                err_msg = "[%s] section is not set" % (self.SECTION_READGROUP)
                raise ValueError(err_msg)
            parsed_bam_tofastq = self.parse_data_readgroup(splited[self.SECTION_READGROUP], bwa_samples, self.SECTION_READGROUP)
            self.readgroup.update(parsed_bam_tofastq["metadata"])
            self.readgroup_src.update(parsed_bam_tofastq["metadata_src"])
