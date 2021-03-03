#! /usr/bin/env python

import gcat_workflow_cloud.abstract_task as abstract_task

class Task(abstract_task.Abstract_task):
    CONF_SECTION = "gatk_mutectcaller_parabricks_compatible"
    TASK_NAME = CONF_SECTION

    def __init__(self, task_dir, sample_conf, param_conf, run_conf):

        super(Task, self).__init__(
            "compat-mutectcaller.sh",
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
                    "--input-recursive REFERENCE_DIR",
                    "--env REFERENCE_FASTA",
                    "--input INPUT_TUMOR_CRAM",
                    "--input INPUT_TUMOR_CRAI",
                    "--input INPUT_NORMAL_CRAM",
                    "--input INPUT_NORMAL_CRAI",
                    "--output-recursive OUTPUT_DIR",
                    "--env TUMOR_SAMPLE",
                    "--env NORMAL_SAMPLE",
                    "--env GATK_JAR",
                    "--env MUTECT_JAVA_OPTION",
                    "--env MUTECT_OPTION",
                    "--env NPROC",
                ]) + "\n"
            )
            for (tumor, normal) in sample_conf.mutect_call:
                normal_bam = ""
                normal_bai = ""
                normal_sample = ""
                if normal != None:
                    normal_bam = "%s/cram/%s/%s.markdup.cram" % (run_conf.output_dir, normal, normal)
                    normal_bai = "%s/cram/%s/%s.markdup.cram.crai" % (run_conf.output_dir, normal, normal)
                    normal_sample = "%s" % (normal)

                hout.write(
                    '\t'.join([
                        param_conf.get(self.CONF_SECTION, "reference_dir"),
                        param_conf.get(self.CONF_SECTION, "reference_file"),
                        "%s/cram/%s/%s.markdup.cram" % (run_conf.output_dir, tumor, tumor),
                        "%s/cram/%s/%s.markdup.cram.crai" % (run_conf.output_dir, tumor, tumor),
                        normal_bam,
                        normal_bai,
                        "%s/mutectcaller/%s" % (run_conf.output_dir, tumor),
                        tumor,
                        normal_sample,
                        param_conf.get(self.CONF_SECTION, "gatk_jar"),
                        param_conf.get(self.CONF_SECTION, "mutect_java_option"),
                        param_conf.get(self.CONF_SECTION, "mutect_option"),
                        param_conf.get(self.CONF_SECTION, "mutect_threads_option"),
                    ]) + "\n"
                )

        return task_file
