#! /usr/bin/env python

import gcat_workflow_cloud.abstract_task as abstract_task

class Task(abstract_task.Abstract_task):
    CONF_SECTION = "genomonsv"
    TASK_NAME = CONF_SECTION

    def __init__(self, task_dir, sample_conf, param_conf, run_conf):

        super(Task, self).__init__(
            "genomonsv.sh",
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
                    "--env TUMOR_SAMPLE",
                    "--env NORMAL_SAMPLE",
                    "--env CONTROL_PANEL",
                    "--input-recursive TUMOR_BAM_DIR",
                    "--env TUMOR_BAM",
                    "--input-recursive NORMAL_BAM_DIR",
                    "--env NORMAL_BAM",
                    "--input-recursive REFERENCE_DIR",
                    "--env REFERENCE_FILE",
                    "--input-recursive MERGED_JUNCTION",
                    "--output-recursive TUMOR_OUTPUT_DIR",
                    "--output-recursive NORMAL_OUTPUT_DIR",
                    "--env GENOMONSV_PARSE_OPTION",
                    "--env GENOMONSV_FILT_OPTION",
                    "--env SV_UTILS_FILT_OPTION",
                ]) + "\n"
            )
            for (tumor, normal, controlpanel) in sample_conf.genomon_sv:
                normal_sample = "None"
                normal_bam_dir = ""
                normal_bam = ""
                normal_output_dir = ""
                if normal != None:
                    normal_sample = normal
                    normal_bam_dir = "%s/cram/%s" % (run_conf.output_dir, normal)
                    normal_bam = "%s.markdup.cram" % (normal)
                    normal_output_dir = "%s/genomonsv/%s" % (run_conf.output_dir, normal),
                
                hout.write(
                    '\t'.join([
                        tumor,
                        normal_sample,
                        "None",
                        "%s/cram/%s" % (run_conf.output_dir, tumor),
                        "%s.markdup.cram" % (tumor),
                        normal_bam_dir,
                        normal_bam,
                        param_conf.get(self.CONF_SECTION, "reference_dir"),
                        param_conf.get(self.CONF_SECTION, "reference_file"),
                        "",
                        "%s/genomonsv/%s" % (run_conf.output_dir, tumor),
                        normal_output_dir,
                        param_conf.get(self.CONF_SECTION, "genomonsv_parse_option"),
                        param_conf.get(self.CONF_SECTION, "genomonsv_filt_option"),
                        param_conf.get(self.CONF_SECTION, "sv_utils_filt_option"),
                    ]) + "\n"
                )

        return task_file
