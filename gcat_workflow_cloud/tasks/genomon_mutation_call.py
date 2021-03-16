#! /usr/bin/env python

import gcat_workflow_cloud.abstract_task as abstract_task

class Task(abstract_task.Abstract_task):
    CONF_SECTION = "genomon_mutation_call"
    TASK_NAME = CONF_SECTION

    def __init__(self, task_dir, sample_conf, param_conf, run_conf):

        super(Task, self).__init__(
            "genomon_mutation_call.sh",
            param_conf.get(self.CONF_SECTION, "image"),
            param_conf.get(self.CONF_SECTION, "resource"),
            run_conf.output_dir + "/logging"
        )
        self.task_file = self.task_file_generation(task_dir, sample_conf, param_conf, run_conf)
        

    def task_file_generation(self, task_dir, sample_conf, param_conf, run_conf):
        task_file = "{}/{}-tasks-{}.tsv".format(task_dir, self.TASK_NAME, run_conf.project_name)
        with open(task_file, 'w') as hout:
            
            hout.write(
                '\t'.join([
                    "--input-recursive INPUT_BAM_DIR1",
                    "--env INPUT_BAM1",
                    "--input-recursive INPUT_BAM_DIR2",
                    "--env INPUT_BAM2",
                    "--input REFERENCE",
                    "--input REFERENCE_IDX",
                    "--input FISHER_INTERVAL_LIST",
                    "--env FISHER_PAIR_OPTION",
                    "--env FISHER_SINGLE_OPTION",
                    "--env FISHER_SAMTOOLS_OPTION",
                    "--env REALIGNMENT_OPTION",
                    "--env INDEL_OPTION",
                    "--env INDEL_SAMTOOLS_OPTION",
                    "--env BREAKPOINT_OPTION",
                    "--output-recursive OUTPUT_DIR",
                    "--env SAMPLE1",
                    "--env SAMPLE2",
                    "--input-recursive ANNOTATION_DB",
                    "--env FILTER_PAIR_OPTION",
                    "--env FILTER_SINGLE_OPTION",
                ]) + "\n"
            )
            for (tumor, normal, controlpanel) in sample_conf.genomon_mutation_call:
                normal_sample = "None"
                normal_bam_dir = ""
                normal_bam = ""
                if normal != None:
                    normal_sample = normal
                    normal_bam_dir = "{root}/{cram}/{sample}".format(root = run_conf.output_dir, sample = normal, cram = run_conf.seq_format)
                    normal_bam = "{sample}.markdup.{cram}".format(sample = normal, cram = run_conf.seq_format)
                
                hout.write(
                    '\t'.join([
                        "{root}/{cram}/{sample}".format(root = run_conf.output_dir, sample = tumor, cram = run_conf.seq_format),
                        "{sample}.markdup.{cram}".format(sample = tumor, cram = run_conf.seq_format),
                        normal_bam_dir,
                        normal_bam,
                        param_conf.get(self.CONF_SECTION, "reference"),
                        param_conf.get(self.CONF_SECTION, "reference_idx"),
                        param_conf.get(self.CONF_SECTION, "fisher_interval_list"),
                        param_conf.get(self.CONF_SECTION, "fisher_pair_option"),
                        param_conf.get(self.CONF_SECTION, "fisher_single_option"),
                        param_conf.get(self.CONF_SECTION, "fisher_samtools_option"),
                        param_conf.get(self.CONF_SECTION, "realignment_option"),
                        param_conf.get(self.CONF_SECTION, "indel_option"),
                        param_conf.get(self.CONF_SECTION, "indel_samtools_option"),
                        param_conf.get(self.CONF_SECTION, "breakpoint_option"),
                        "%s/genomon_mutation_call/%s" % (run_conf.output_dir, tumor),
                        tumor,
                        normal_sample,
                        param_conf.get(self.CONF_SECTION, "annotation_db"),
                        param_conf.get(self.CONF_SECTION, "filter_pair_option"),
                        param_conf.get(self.CONF_SECTION, "filter_single_option"),
                    ]) + "\n"
                )

        return task_file
