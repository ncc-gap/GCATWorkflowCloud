#! /usr/bin/env python

import gcat_workflow_cloud.abstract_task as abstract_task

class Task(abstract_task.Abstract_task):
    CONF_SECTION = "manta"
    TASK_NAME = CONF_SECTION

    def __init__(self, task_dir, sample_conf, param_conf, run_conf):

        if run_conf.analysis_type == "germline":
            script_name = "manta-germline.sh"

        elif run_conf.analysis_type == "somatic":
            script_name = "manta-somatic.sh"

        super(Task, self).__init__(
            script_name,
            param_conf.get(self.CONF_SECTION, "image"),
            param_conf.get(self.CONF_SECTION, "resource"),
            run_conf.output_dir + "/logging"
        )
        
        self.task_file = self.task_file_generation(task_dir, sample_conf, param_conf, run_conf)

    def _germline(self, task_dir, sample_conf, param_conf, run_conf):

        task_file = "{}/{}-tasks-{}.tsv".format(task_dir, self.TASK_NAME, run_conf.project_name)
        with open(task_file, 'w') as hout:
            
            hout.write(
                '\t'.join([
                    "--input-recursive NORMAL_BAM_DIR",
                    "--env NORMAL_BAM_FILE",
                    "--output-recursive OUTPUT_DIR",
                    "--input REFERENCE",
                    "--input REFERENCE_IDX",
                    "--env CONFIG_MANTA_OPTION",
                    "--env WORKFLOW_OPTION",
                ]) + "\n"
            )
            for sample in sample_conf.manta:
                hout.write(
                    '\t'.join([
                        "{root}/{cram}/{sample}".format(root = run_conf.output_dir, sample = sample, cram = run_conf.seq_format),
                        "{sample}.markdup.{cram}".format(sample = sample, cram = run_conf.seq_format),
                        "%s/manta/%s" % (run_conf.output_dir, sample),
                        param_conf.get(self.CONF_SECTION, "reference"),
                        param_conf.get(self.CONF_SECTION, "reference_idx"),
                        param_conf.get(self.CONF_SECTION, "manta_config_option"),
                        param_conf.get(self.CONF_SECTION, "manta_workflow_option") + " " + param_conf.get(self.CONF_SECTION, "manta_workflow_threads_option"),
                    ]) + "\n"
                )
        return task_file

    def _somatic(self, task_dir, sample_conf, param_conf, run_conf):

        task_file = "{}/{}-tasks-{}.tsv".format(task_dir, self.TASK_NAME, run_conf.project_name)
        with open(task_file, 'w') as hout:
            
            hout.write(
                '\t'.join([
                    "--input-recursive NORMAL_BAM_DIR",
                    "--env NORMAL_BAM_FILE",
                    "--input-recursive TUMOR_BAM_DIR",
                    "--env TUMOR_BAM_FILE",
                    "--output-recursive OUTPUT_DIR",
                    "--input REFERENCE",
                    "--input REFERENCE_IDX",
                    "--env CONFIG_MANTA_OPTION",
                    "--env WORKFLOW_OPTION",
                ]) + "\n"
            )
            for (tumor, normal) in sample_conf.manta:
                normal_bam_dir = ""
                normal_bam_file = ""
                if normal != None:
                    normal_bam_dir = "{root}/{cram}/{sample}".format(root = run_conf.output_dir, sample = normal, cram = run_conf.seq_format)
                    normal_bam_file = "{sample}.markdup.{cram}".format(sample = normal, cram = run_conf.seq_format)
                    
                hout.write(
                    '\t'.join([
                        normal_bam_dir,
                        normal_bam_file,
                        "{root}/{cram}/{sample}".format(root = run_conf.output_dir, sample = tumor, cram = run_conf.seq_format),
                        "{sample}.markdup.{cram}".format(sample = tumor, cram = run_conf.seq_format),
                        "%s/manta/%s" % (run_conf.output_dir, tumor),
                        param_conf.get(self.CONF_SECTION, "reference"),
                        param_conf.get(self.CONF_SECTION, "reference_idx"),
                        param_conf.get(self.CONF_SECTION, "manta_config_option"),
                        param_conf.get(self.CONF_SECTION, "manta_workflow_option") + " " + param_conf.get(self.CONF_SECTION, "manta_workflow_threads_option"),
                    ]) + "\n"
                )
        return task_file

    def task_file_generation(self, task_dir, sample_conf, param_conf, run_conf):
        if run_conf.analysis_type == "germline":
            return self._germline(task_dir, sample_conf, param_conf, run_conf)

        elif run_conf.analysis_type == "somatic":
            return self._somatic(task_dir, sample_conf, param_conf, run_conf)

        return None
