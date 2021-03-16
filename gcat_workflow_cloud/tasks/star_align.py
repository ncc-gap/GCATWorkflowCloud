#! /usr/bin/env python

import gcat_workflow_cloud.abstract_task as abstract_task

class Task(abstract_task.Abstract_task):
    CONF_SECTION = "star_alignment"
    TASK_NAME = CONF_SECTION

    def __init__(self, task_dir, sample_conf, param_conf, run_conf):

        super(Task, self).__init__(
            "star-align.sh",
            param_conf.get(self.CONF_SECTION, "image"),
            param_conf.get(self.CONF_SECTION, "resource"),
            run_conf.output_dir + "/logging"
        )
        
        self.task_file = self.task_file_generation(task_dir, sample_conf, param_conf, run_conf)

    def task_file_generation(self, task_dir, sample_conf, param_conf, run_conf):

        task_file = "{}/{}-tasks-{}.tsv".format(task_dir, self.TASK_NAME, run_conf.project_name)

        max_len_fq1 = 0
        max_len_fq2 = 0
        for sample in sample_conf.fastq:
            for i,fq in enumerate(sample_conf.fastq[sample]):
                if i == 0 and len(fq) > max_len_fq1:
                    max_len_fq1 = len(fq)
                if i == 1 and len(fq) > max_len_fq2:
                    max_len_fq2 = len(fq)

        header_fq1 = []
        for i in range(max_len_fq1):
            header_fq1.append("--input FASTQ1_%d" % (i+1))

        header_fq2 = []
        for i in range(max_len_fq2):
            header_fq2.append("--input FASTQ2_%d" % (i+1))

        with open(task_file, 'w') as hout:
            
            hout.write(
                '\t'.join(header_fq1 + header_fq2 + [
                    "--output-recursive OUTPUT_DIR",
                    "--input-recursive STAR_REFERENCE",
                    "--input REFERENCE",
                    "--input REFERENCE_INDEX",
                    "--env STAR_OPTION",
                    "--env SAMPLE",
                    "--env FQ_TYPE",
                    "--env SEQ_FORMAT",
                ]) + "\n"
            )
            for sample in sample_conf.fastq:
                fq1 = (sample_conf.fastq[sample][0] + [""] * max_len_fq1)[0:max_len_fq1]
                fq2_empty = [""] * max_len_fq2
                fq_type = "single"
                if len(sample_conf.fastq[sample]) > 1:
                    fq2 = (sample_conf.fastq[sample][1] + fq2_empty)[0:max_len_fq2]
                    fq_type = "pair"
                hout.write(
                    '\t'.join(fq1 + fq2 + [
                        "{root}/{cram}/{sample}".format(root = run_conf.output_dir, sample = sample, cram = run_conf.seq_format),
                        param_conf.get(self.CONF_SECTION, "star_genome"),
                        param_conf.get(self.CONF_SECTION, "reference"),
                        param_conf.get(self.CONF_SECTION, "reference_index"),
                        param_conf.get(self.CONF_SECTION, "star_option"),
                        sample,
                        fq_type,
                        run_conf.seq_format,
                    ]) + "\n"
                )
        return task_file
