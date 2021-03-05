#! /usr/bin/env python

import gcat_workflow_cloud.abstract_task as abstract_task

TAGS = {
    "chr1":        {"ploidy": "2", "interval": "interval_file_autosome_chr1", "label": "autosome.chr1"},
    "chr2":        {"ploidy": "2", "interval": "interval_file_autosome_chr2", "label": "autosome.chr2"},
    "chr3":        {"ploidy": "2", "interval": "interval_file_autosome_chr3", "label": "autosome.chr3"},
    "chr4":        {"ploidy": "2", "interval": "interval_file_autosome_chr4", "label": "autosome.chr4"},
    "chr5":        {"ploidy": "2", "interval": "interval_file_autosome_chr5", "label": "autosome.chr5"},
    "chr6":        {"ploidy": "2", "interval": "interval_file_autosome_chr6", "label": "autosome.chr6"},
    "chr7":        {"ploidy": "2", "interval": "interval_file_autosome_chr7", "label": "autosome.chr7"},
    "chr8":        {"ploidy": "2", "interval": "interval_file_autosome_chr8", "label": "autosome.chr8"},
    "chr9":        {"ploidy": "2", "interval": "interval_file_autosome_chr9", "label": "autosome.chr9"},
    "chr10":       {"ploidy": "2", "interval": "interval_file_autosome_chr10", "label": "autosome.chr10"},
    "chr11":       {"ploidy": "2", "interval": "interval_file_autosome_chr11", "label": "autosome.chr11"},
    "chr12":       {"ploidy": "2", "interval": "interval_file_autosome_chr12", "label": "autosome.chr12"},
    "chr13":       {"ploidy": "2", "interval": "interval_file_autosome_chr13", "label": "autosome.chr13"},
    "chr14":       {"ploidy": "2", "interval": "interval_file_autosome_chr14", "label": "autosome.chr14"},
    "chr15":       {"ploidy": "2", "interval": "interval_file_autosome_chr15", "label": "autosome.chr15"},
    "chr16":       {"ploidy": "2", "interval": "interval_file_autosome_chr16", "label": "autosome.chr16"},
    "chr17":       {"ploidy": "2", "interval": "interval_file_autosome_chr17", "label": "autosome.chr17"},
    "chr18":       {"ploidy": "2", "interval": "interval_file_autosome_chr18", "label": "autosome.chr18"},
    "chr19":       {"ploidy": "2", "interval": "interval_file_autosome_chr19", "label": "autosome.chr19"},
    "chr20":       {"ploidy": "2", "interval": "interval_file_autosome_chr20", "label": "autosome.chr20"},
    "chr21":       {"ploidy": "2", "interval": "interval_file_autosome_chr21", "label": "autosome.chr21"},
    "chr22":       {"ploidy": "2", "interval": "interval_file_autosome_chr22", "label": "autosome.chr22"},
    "chrx_female": {"ploidy": "2", "interval": "interval_file_chrx", "label": "chrX.female"},
    "chrx_male":   {"ploidy": "1", "interval": "interval_file_chrx", "label": "chrX.male"},
    "chry_male":   {"ploidy": "1", "interval": "interval_file_chry", "label": "chrY.male"},
    "par":         {"ploidy": "2", "interval": "interval_file_par", "label": "PAR"},
}

class Task(abstract_task.Abstract_task):
    CONF_SECTION = "gatk_haplotypecaller_parabricks_compatible"
    TASK_NAME = CONF_SECTION

    def __init__(self, task_dir, sample_conf, param_conf, run_conf, key):

        super(Task, self).__init__(
            "compat-haplotypecaller.sh",
            param_conf.get(self.CONF_SECTION, "image"),
            param_conf.get(self.CONF_SECTION, "resource"),
            run_conf.output_dir + "/logging"
        )
        self.TASK_NAME += "." + key
        self.task_file = self.task_file_generation(task_dir, sample_conf, param_conf, run_conf, key)
        

    def task_file_generation(self, task_dir, sample_conf, param_conf, run_conf, key):
        task_file = "{}/{}-tasks-{}.tsv".format(task_dir, self.TASK_NAME, run_conf.project_name)
        with open(task_file, 'w') as hout:
            
            hout.write(
                '\t'.join([
                    "--input REFERENCE",
                    "--input REFERENCE_IDX",
                    "--input REFERENCE_DICT",
                    "--input INTERVAL",
                    "--input INPUT_CRAM",
                    "--input INPUT_CRAI",
                    "--output-recursive OUTPUT_DIR",
                    "--env SAMPLE",
                    "--env GATK_JAR",
                    "--env HAPLOTYPE_JAVA_OPTION",
                    "--env HAPLOTYPE_OPTION",
                    "--env NPROC",
                    "--env PLOIDY",
                    "--env TAG",
                ]) + "\n"
            )
            for sample in sample_conf.__dict__["haplotype_call_" + key]:
                
                hout.write(
                    '\t'.join([
                        param_conf.get(self.CONF_SECTION, "reference"),
                        param_conf.get(self.CONF_SECTION, "reference_idx"),
                        param_conf.get(self.CONF_SECTION, "reference_dict"),
                        param_conf.get(self.CONF_SECTION, TAGS[key]["interval"]),
                        "%s/cram/%s/%s.markdup.cram" % (run_conf.output_dir, sample, sample),
                        "%s/cram/%s/%s.markdup.cram.crai" % (run_conf.output_dir, sample, sample),
                        "%s/haplotypecaller/%s" % (run_conf.output_dir, sample),
                        "%s" % (sample),
                        param_conf.get(self.CONF_SECTION, "gatk_jar"),
                        param_conf.get(self.CONF_SECTION, "haplotype_java_option"),
                        param_conf.get(self.CONF_SECTION, "haplotype_option"),
                        param_conf.get(self.CONF_SECTION, "haplotype_threads_option"),
                        TAGS[key]["ploidy"],
                        TAGS[key]["label"],
                    ]) + "\n"
                )

        return task_file
