# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:18:34 2016

@author: okada

"""

import unittest
import subprocess
import time

class ConfigureTest(unittest.TestCase):
    
    # init class
    @classmethod
    def setUpClass(self):
        pass

    # terminated class
    @classmethod
    def tearDownClass(self):
        pass

    # init method
    def setUp(self):
        pass

    # terminated method
    def tearDown(self):
        pass
    
    def test1_01_version(self):
        subprocess.check_call(['gcat_workflow_cloud', '--version'])
    
    def test2_01_configure(self):
        options = [
            "germline",
            "./example_conf/germline_sample.csv",
            "s3://travisci-work",
            "./example_conf/germline_gcat_ecsub.cfg",
            "--dryrun",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test2_02_configure(self):
        options = [
            "germline",
            "./example_conf/germline_sample.csv",
            "s3://travisci-work",
            "./example_conf/germline_gcat_ecsub.cfg",
            "--dryrun",
            "--skip-price",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test2_03_configure(self):
        options = [
            "germline",
            "./example_conf/germline_sample.csv",
            "s3://travisci-work",
            "./example_conf/germline_gcat_ecsub.cfg",
            "--dryrun",
            "--use-bam",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test2_04_configure(self):
        options = [
            "germline",
            "./tests/germline_sample_allinone_cram.csv",
            "s3://travisci-work",
            "./example_conf/germline_gcat_ecsub.cfg",
            "--dryrun",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test2_05_configure(self):
        options = [
            "germline",
            "./tests/germline_sample_allinone_bam.csv",
            "s3://travisci-work",
            "./example_conf/germline_gcat_ecsub.cfg",
            "--dryrun",
            "--use-bam",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test2_06_configure(self):
        options = [
            "germline",
            "./tests/germline_sample_fastq.csv",
            "s3://travisci-work",
            "./example_conf/germline_gcat_ecsub.cfg",
            "--dryrun",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test2_07_configure(self):
        options = [
            "germline",
            "./tests/germline_sample_fastq.csv",
            "s3://travisci-work",
            "./example_conf/germline_gcat_ecsub.cfg",
            "--dryrun",
            "--use-bam",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test2_08_configure(self):
        options = [
            "germline",
            "./tests/germline_sample_bamimport_cram.csv",
            "s3://travisci-work",
            "./example_conf/germline_gcat_ecsub.cfg",
            "--dryrun",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test2_09_configure(self):
        options = [
            "germline",
            "./tests/germline_sample_bamimport_bam.csv",
            "s3://travisci-work",
            "./example_conf/germline_gcat_ecsub.cfg",
            "--dryrun",
            "--use-bam",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test2_10_configure(self):
        options = [
            "germline",
            "./tests/germline_sample_bamtofastq.csv",
            "s3://travisci-work",
            "./example_conf/germline_gcat_ecsub.cfg",
            "--dryrun",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test2_11_configure(self):
        options = [
            "germline",
            "./tests/germline_sample_bamtofastq.csv",
            "s3://travisci-work",
            "./example_conf/germline_gcat_ecsub.cfg",
            "--dryrun",
            "--use-bam",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test3_01_configure(self):
        options = [
            "rna",
            "./example_conf/rna_sample.csv",
            "s3://travisci-work",
            "./example_conf/rna_gcat_ecsub.cfg",
            "--dryrun",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)
        time.sleep(2)

    def test3_02_configure(self):
        options = [
            "rna",
            "./example_conf/rna_sample.csv",
            "s3://travisci-work",
            "./example_conf/rna_gcat_ecsub.cfg",
            "--dryrun",
            "--skip-price",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)
        time.sleep(2)

    def test3_03_configure(self):
        options = [
            "rna",
            "./example_conf/rna_sample.csv",
            "s3://travisci-work",
            "./example_conf/rna_gcat_ecsub.cfg",
            "--dryrun",
            "--use-bam",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)
        time.sleep(2)
        
    def test4_01_configure(self):
        options = [
            "somatic",
            "./example_conf/somatic_sample.csv",
            "s3://travisci-work",
            "./example_conf/somatic_gcat_ecsub.cfg",
            "--dryrun",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test4_02_configure(self):
        options = [
            "somatic",
            "./example_conf/somatic_sample.csv",
            "s3://travisci-work",
            "./example_conf/somatic_gcat_ecsub.cfg",
            "--dryrun",
            "--skip-price",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)
        
    def test4_03_configure(self):
        options = [
            "somatic",
            "./example_conf/somatic_sample.csv",
            "s3://travisci-work",
            "./example_conf/somatic_gcat_ecsub.cfg",
            "--dryrun",
            "--use-bam",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test4_04_configure(self):
        options = [
            "somatic",
            "./tests/somatic_sample_allinone_cram.csv",
            "s3://travisci-work",
            "./example_conf/somatic_gcat_ecsub.cfg",
            "--dryrun",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test4_05_configure(self):
        options = [
            "somatic",
            "./tests/somatic_sample_allinone_bam.csv",
            "s3://travisci-work",
            "./example_conf/somatic_gcat_ecsub.cfg",
            "--dryrun",
            "--use-bam",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test4_06_configure(self):
        options = [
            "somatic",
            "./tests/somatic_sample_fastq.csv",
            "s3://travisci-work",
            "./example_conf/somatic_gcat_ecsub.cfg",
            "--dryrun",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test4_07_configure(self):
        options = [
            "somatic",
            "./tests/somatic_sample_fastq.csv",
            "s3://travisci-work",
            "./example_conf/somatic_gcat_ecsub.cfg",
            "--dryrun",
            "--use-bam",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test4_08_configure(self):
        options = [
            "somatic",
            "./tests/somatic_sample_bamimport_cram.csv",
            "s3://travisci-work",
            "./example_conf/somatic_gcat_ecsub.cfg",
            "--dryrun",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test4_09_configure(self):
        options = [
            "somatic",
            "./tests/somatic_sample_bamimport_bam.csv",
            "s3://travisci-work",
            "./example_conf/somatic_gcat_ecsub.cfg",
            "--dryrun",
            "--use-bam",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test4_10_configure(self):
        options = [
            "somatic",
            "./tests/somatic_sample_bamtofastq.csv",
            "s3://travisci-work",
            "./example_conf/somatic_gcat_ecsub.cfg",
            "--dryrun",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

    def test4_11_configure(self):
        options = [
            "somatic",
            "./tests/somatic_sample_bamtofastq.csv",
            "s3://travisci-work",
            "./example_conf/somatic_gcat_ecsub.cfg",
            "--dryrun",
            "--use-bam",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

if __name__ == '__main__':
    unittest.main()
