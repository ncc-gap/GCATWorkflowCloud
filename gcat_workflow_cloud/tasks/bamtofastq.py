#! /usr/bin/env python

import gcat_workflow_cloud.abstract_task as abstract_task

class Task(abstract_task.Abstract_task):
    CONF_SECTION = "bam_tofastq"
    TASK_NAME = CONF_SECTION

    def __init__(self, task_dir, sample_conf, param_conf, run_conf):

        super(Task, self).__init__(
            "bamtofastq.sh",
            param_conf.get(self.CONF_SECTION, "image"),
            param_conf.get(self.CONF_SECTION, "resource"),
            run_conf.output_dir + "/logging"
        )
        
        self.task_file = self.task_file_generation(task_dir, sample_conf, param_conf, run_conf)

    def task_file_generation(self, task_dir, sample_conf, param_conf, run_conf):

        task_file = "{}/{}-tasks-{}.tsv".format(task_dir, self.TASK_NAME, run_conf.project_name)

        max_len_bam = 0
        for sample in sample_conf.bam_tofastq:
            for i,path in enumerate(sample_conf.bam_tofastq[sample]):
                if i == 0 and len(path) > max_len_bam:
                    max_len_bam = len(path)


        header_bam = []
        for i in range(max_len_bam):
            header_bam.append("--input INPUT_BAM_%d" % (i+1))

        with open(task_file, 'w') as hout:
            
            hout.write(
                '\t'.join(header_bam + [
                    "--output-recursive OUTPUT_DIR",
                    "--env SAMPLE_MAX_INDEX",
                    "--env PARAM",
                ]) + "\n"
            )
            for sample in sample_conf.bam_tofastq:
                bam = ([sample_conf.bam_tofastq[sample]] + [""] * max_len_bam)[0:max_len_bam]

                hout.write(
                    '\t'.join(bam + [
                        "%s/fastq/%s" % (run_conf.output_dir, sample),
                        str(max_len_bam),
                        param_conf.get(self.CONF_SECTION, "option"),
                    ]) + "\n"
                )
        return task_file
