# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:18:34 2016

@author: okada

"""

import unittest
import subprocess

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
    
    def test2_01_configure_drmaa_nogpu(self):
        options = [
            "germline",
            "./tests/germline_sample.csv",
            "s3://travisci-work",
            "./tests/germline_gcat_ecsub.cfg",
            "--dryrun",
        ]
        subprocess.check_call(['gcat_workflow_cloud'] + options)

if __name__ == '__main__':
    unittest.main()
