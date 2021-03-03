#! /usr/bin/env python

import gcat_workflow_cloud.abstract_task as abstract_task

class Task(abstract_task.Abstract_task):
    CONF_SECTION = "genomon_mutation_call"
    TASK_NAME = CONF_SECTION

    def __init__(self, task_dir, sample_conf, param_conf, run_conf):

        super(Task, self).__init__(
            "genomon_mutation_call.sh",
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
                ]) + "\n"
            )
"""
${ACTIVE_EXAC_FLAG}
${ACTIVE_HGVD_2016_FLAG}
${ANNOTATION_DB}
${BLAT} 
${BREAKPOINT_OPTION}
${CONTROL_BAM_LIST}
${FILTER_PAIR_OPTION}
${FILTER_SINGLE_OPTION}
${FISHER_INTERVAL_LIST}
${FISHER_PAIR_OPTION}
${FISHER_PAIR_SAMTOOLS}
${FISHER_SINGLE_OPTION}
${FISHER_SINGLE_SAMTOOLS}
${HOTSPOT_DB}
${HOTSPOT_OPTION}
${HOTSPOT_SAMTOOLS}
${INDEL_OPTION}
${INDEL_SAMTOOLS}
${INPUT_BAM1}
${INPUT_BAM2}
${INPUT_DIR1}
${INPUT_DIR2}
${OUTPUT_DIR}
${PATH}
${REALIGNMENT_OPTION}
${REFERENCE}
${SAMPLE1}
${SAMPLE2}
${SAMTOOLS} 
"""
            for (tumor, normal, controlpanel) in sample_conf.genomon_mutation_call:
                
                hout.write(
                    '\t'.join([
                        param_conf.get(self.CONF_SECTION, "reference_dir"),
                        param_conf.get(self.CONF_SECTION, "reference_file"),
                        "%s/cram/%s/%s.markdup.cram" % (run_conf.output_dir, tumor, tumor),
                        "%s/cram/%s/%s.markdup.cram.crai" % (run_conf.output_dir, tumor, tumor),
                        "%s/cram/%s/%s.markdup.cram" % (run_conf.output_dir, normal, normal),
                        "%s/cram/%s/%s.markdup.cram.crai" % (run_conf.output_dir, normal, normal),
                        "%s/genomon_mutation_call/%s" % (run_conf.output_dir, tumor),
                        "%s" % (tumor),
                        "%s" % (normal),
                    ]) + "\n"
                )

        return task_file
