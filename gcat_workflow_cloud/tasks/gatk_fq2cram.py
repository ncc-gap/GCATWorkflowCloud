#! /usr/bin/env python

import gcat_workflow_cloud.abstract_task as abstract_task

class Task(abstract_task.Abstract_task):
    CONF_SECTION = "gatk_bwa_alignment_parabricks_compatible"
    TASK_NAME = CONF_SECTION

    def __init__(self, output_dir, task_dir, sample_conf, param_conf, run_conf):
        
        super(Task, self).__init__(
            "compat-fq2cram.sh",
            param_conf.get(self.CONF_SECTION, "image"),
            param_conf.get(self.CONF_SECTION, "resource"),
            output_dir + "/logging"
        )

        self.task_file = self.task_file_generation(output_dir, task_dir, sample_conf, param_conf, run_conf)

    def task_file_generation(self, output_dir, task_dir, sample_conf, param_conf, run_conf):

        task_file = "{}/{}-tasks-{}.tsv".format(task_dir, self.TASK_NAME, run_conf.project_name)
        with open(task_file, 'w') as hout:
            hout.write(
                '\t'.join([
                    "--input-recursive REFERENCE_DIR",
                    "--env REFERENCE_FASTA",
                    "--env SAMPLE_NAME",
                    "--input INPUT_FASTQ_1",
                    "--input INPUT_FASTQ_2",
                    "--output OUTPUT_CRAM",
                    "--output OUTPUT_CRAI",
                    "--output OUTPUT_MARKDUP_METRICS",
                ]) + "\n"
            )
            for sample in sample_conf.fastq:
                        
                if len(sample_conf.fastq[sample][0]) == 1:
                    fastq1 = sample_conf.fastq[sample][0][0]
                else:
                    fastq1 = "'<cat %s'" % (" ".join(sample_conf.fastq[sample][0]))

                if len(sample_conf.fastq[sample][1]) == 1:
                    fastq2 = sample_conf.fastq[sample][1][0]
                else:
                    fastq2 = "'<cat %s'" % (" ".join(sample_conf.fastq[sample][1]))

                hout.write(
                    '\t'.join([
                        param_conf.get(self.CONF_SECTION, "reference_dir"),
                        param_conf.get(self.CONF_SECTION, "reference_file"),
                        sample,
                        fastq1,
                        fastq2,
                        "%s/cram/%s/%s.markdup.cram" % (output_dir, sample, sample),
                        "%s/cram/%s/%s.markdup.cram.crai" % (output_dir, sample, sample),
                        "%s/cram/%s/%s.markdup.metrics" % (output_dir, sample, sample),
                    ]) + "\n"
                )
        return task_file
