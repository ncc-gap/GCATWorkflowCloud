#! /usr/bin/env python

import gcat_workflow_cloud.abstract_task as abstract_task

class Task(abstract_task.Abstract_task):
    CONF_SECTION = "gridss"
    TASK_NAME = CONF_SECTION

    def __init__(self, task_dir, sample_conf, param_conf, run_conf):

        if run_conf.analysis_type == "germline":
            script_name = "gridss-germline.sh"

        elif run_conf.analysis_type == "somatic":
            script_name = "gridss-somatic.sh",

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
                    "--output-recursive OUTPUT_DIR",
                    "--env VCF",
                    "--env ASSEMBLE",
                    "--input-recursive REFERENCE_DIR",
                    "--env REFERENCE_FILE",
                    "--input INPUT_CRAM",
                    "--input INPUT_CRAI",
                    "--env GRIDSS_JAR",
                ]) + "\n"
            )
            for sample in sample_conf.gridss:
                hout.write(
                    '\t'.join([
                        "%s/gridss/%s" % (run_conf.output_dir, sample),
                        "%s_gridss-result.vcf" % (sample),
                        "%s_gridss-assembly.bam" % (sample),
                        param_conf.get(self.CONF_SECTION, "reference_dir"),
                        param_conf.get(self.CONF_SECTION, "reference_file"),
                        "%s/cram/%s/%s.markdup.cram" % (run_conf.output_dir, sample, sample),
                        "%s/cram/%s/%s.markdup.cram.crai" % (run_conf.output_dir, sample, sample),
                        param_conf.get(self.CONF_SECTION, "gridss_jar"),
                    ]) + "\n"
                )

        return task_file

    def _somatic(self, task_dir, sample_conf, param_conf, run_conf):

        task_file = "{}/{}-tasks-{}.tsv".format(task_dir, self.TASK_NAME, run_conf.project_name)
        with open(task_file, 'w') as hout:
            
            hout.write(
                '\t'.join([
                    "--output-recursive OUTPUT_DIR",
                    "--env VCF",
                    "--env VCF_SOMATIC",
                    "--env ASSEMBLE",
                    "--input-recursive REFERENCE_DIR",
                    "--env REFERENCE_FILE",
                    "--input NORMAL_CRAM",
                    "--input NORMAL_CRAI",
                    "--input TUMOR_CRAM",
                    "--input TUMOR_CRAI",
                    "--env GRIDSS_JAR",
                ]) + "\n"
            )

            for (tumor, normal) in sample_conf.gridss:
                normal_bam = ""
                normal_bai = ""
                if normal != None:
                    normal_bam = "%s/cram/%s/%s.markdup.cram" % (run_conf.output_dir, normal, normal)
                    normal_bai = "%s/cram/%s/%s.markdup.cram.crai" % (run_conf.output_dir, normal, normal)

                hout.write(
                    '\t'.join([
                        "%s/gridss/%s" % (run_conf.output_dir, tumor),
                        "%s_gridss-result.vcf" % (tumor),
                        "%s_gridss-result.somatic.vcf" % (tumor),
                        "%s_gridss-assembly.bam" % (tumor),
                        param_conf.get(self.CONF_SECTION, "reference_dir"),
                        param_conf.get(self.CONF_SECTION, "reference_file"),
                        normal_bam,
                        normal_bai,
                        "%s/cram/%s/%s.markdup.cram" % (run_conf.output_dir, tumor, tumor),
                        "%s/cram/%s/%s.markdup.cram.crai" % (run_conf.output_dir, tumor, tumor),
                        param_conf.get(self.CONF_SECTION, "gridss_jar"),
                    ]) + "\n"
                )

        return task_file

    def task_file_generation(self, task_dir, sample_conf, param_conf, run_conf):
        if run_conf.analysis_type == "germline":
            return self._germline(task_dir, sample_conf, param_conf, run_conf)

        elif run_conf.analysis_type == "somatic":
            return self._somatic(task_dir, sample_conf, param_conf, run_conf)

        return None
