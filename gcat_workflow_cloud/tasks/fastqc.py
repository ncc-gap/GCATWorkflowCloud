#! /usr/bin/env python

import gcat_workflow_cloud.abstract_task as abstract_task

class Task(abstract_task.Abstract_task):
    CONF_SECTION = "fastqc"
    TASK_NAME = CONF_SECTION

    def __init__(self, output_dir, task_dir, sample_conf, param_conf, run_conf):

        super(Task, self).__init__(
            "fastqc.sh",
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
                    "--input INPUT_FASTQ_R1",
                    "--input INPUT_FASTQ_R2",
                    "--output-recursive OUTPUT_DIR",
                    "--env FASTQC_PARAMS",
                ]) + "\n"
            )
            for sample in sample_conf.fastqc:
                if not sample in sample_conf.fastq:
                    err_msg = "[fastqc] section, %s is not defined in [fastq] section" % (sample)
                    raise ValueError(err_msg)

                hout.write(
                    '\t'.join([
                        ' '.join(sample_conf.fastqc[sample][0]),
                        ' '.join(sample_conf.fastqc[sample][1]),
                        "%s/fastqc/%s" % (output_dir, sample),
                        param_conf.get(self.CONF_SECTION, "fastqc_option"),
                    ]) + "\n"
                )

        return task_file
