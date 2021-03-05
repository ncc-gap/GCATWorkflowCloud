#! /usr/bin/env python

import gcat_workflow_cloud.abstract_task as abstract_task
import os

class Task(abstract_task.Abstract_task):
    CONF_SECTION = "gatk_bwa_alignment_parabricks_compatible"
    TASK_NAME = CONF_SECTION

    def __init__(self, task_dir, sample_conf, param_conf, run_conf):
        
        super(Task, self).__init__(
            "compat-fq2cram.sh",
            param_conf.get(self.CONF_SECTION, "image"),
            param_conf.get(self.CONF_SECTION, "resource"),
            run_conf.output_dir + "/logging"
        )

        self.task_file = self.task_file_generation(task_dir, sample_conf, param_conf, run_conf)

    def task_file_generation(self, task_dir, sample_conf, param_conf, run_conf):

        task_file = "{}/{}-tasks-{}.tsv".format(task_dir, self.TASK_NAME, run_conf.project_name)
        
        input_num = 0
        for sample in sample_conf.fastq:
            if len(sample_conf.fastq[sample][0]) != len(sample_conf.fastq[sample][1]):
                raise ValueError("The number of files does not match between R1 and R2. %s" % sample)
            if input_num < len(sample_conf.fastq[sample][0]):
                input_num = len(sample_conf.fastq[sample][0])

        input_fq1_header = []
        input_fq2_header = []
        for i in range(input_num):
            input_fq1_header.append("--input INPUT_FASTQ_1_%d" % (i))
            input_fq2_header.append("--input INPUT_FASTQ_2_%d" % (i))

        with open(task_file, 'w') as hout:
            hout.write(
                '\t'.join([
                    "--input-recursive REFERENCE_DIR",
                    "--env REFERENCE_FASTA",
                    "--env SAMPLE_NAME",
                    "--env GATK_JAR",
                    "--output OUTPUT_CRAM",
                    "--output OUTPUT_CRAI",
                    "--output OUTPUT_MARKDUP_METRICS",
                    "\t".join(input_fq1_header),
                    "\t".join(input_fq2_header),
                    "--env ARRAY_RG",
                    "--env SAMPLE_MAX_INDEX",
                ]) + "\n"
            )
            for sample in sample_conf.fastq:
                input_fq1 = [""] * input_num
                input_fq2 = [""] * input_num
                array_rg = []
                readgroups = open(sample_conf.readgroup_local[sample]).readlines()
                for i, fq1 in enumerate(sample_conf.fastq[sample][0]):
                    #print((fq1, i))
                    input_fq1[i] = fq1
                    input_fq2[i] = sample_conf.fastq[sample][1][i]
                    array_rg.append(readgroups[i].rstrip())

                hout.write(
                    '\t'.join([
                        param_conf.get(self.CONF_SECTION, "reference_dir"),
                        param_conf.get(self.CONF_SECTION, "reference_file"),
                        sample,
                        param_conf.get(self.CONF_SECTION, "gatk_jar"),
                        "%s/cram/%s/%s.markdup.cram" % (run_conf.output_dir, sample, sample),
                        "%s/cram/%s/%s.markdup.cram.crai" % (run_conf.output_dir, sample, sample),
                        "%s/cram/%s/%s.markdup.metrics" % (run_conf.output_dir, sample, sample),
                        "\t".join(input_fq1),
                        "\t".join(input_fq2),
                        '%s' % (" ".join(array_rg)),
                        str(len(sample_conf.fastq[sample][0]) - 1),
                    ]) + "\n"
                )
        return task_file
