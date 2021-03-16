#! /usr/bin/env python

import gcat_workflow_cloud.abstract_task as abstract_task

class Task(abstract_task.Abstract_task):
    CONF_SECTION = "intron_retention"
    TASK_NAME = CONF_SECTION

    def __init__(self, task_dir, sample_conf, param_conf, run_conf):

        super(Task, self).__init__(
            "intron-retention.sh",
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
                    "--input INPUT_CRAM",
                    "--output OUTPUT_FILE",
                    "--output OUTPUT_INDEX",
                    "--input REFERENCE",
                    "--input REFERENCE_INDEX",
                ]) + "\n"
            )
            for sample in sample_conf.ir_count:
                hout.write(
                    '\t'.join([
                        "{root}/{cram}/{sample}/{sample}.Aligned.sortedByCoord.out.{cram}".format(root = run_conf.output_dir, sample = sample, cram = run_conf.seq_format),
                        "%s/ir_count/%s/%s.ir_simple_count.txt.gz" % (run_conf.output_dir, sample, sample),
                        "%s/ir_count/%s/%s.ir_simple_count.txt.gz.tbi" % (run_conf.output_dir, sample, sample),
                        param_conf.get(self.CONF_SECTION, "reference"),
                        param_conf.get(self.CONF_SECTION, "reference_index"),
                    ]) + "\n"
                )
        return task_file
