#! /usr/bin/env python

import os
import datetime

date_format = "{year:0>4d}-{month:0>2d}-{day:0>2d} {hour:0>2d}:{min:0>2d}:{second:0>2d}"
timestamp_format = "{year:0>4d}{month:0>2d}{day:0>2d}-{hour:0>2d}{min:0>2d}{second:0>2d}"

import gcat_workflow_cloud.__init__ as ini

class Run_conf(object):
    """
    class for job related parameters
    """

    def __init__(self, sample_conf_file = None, 
                       param_conf_file = None,
                       analysis_type = None,
                       output_dir = None):

        self.sample_conf_file = sample_conf_file
        self.sample_conf_name = os.path.splitext(os.path.basename(sample_conf_file))[0]
        self.param_conf_file = param_conf_file
        self.analysis_type = analysis_type
        
        now = datetime.datetime.now()
        # format "YYYY-mm-dd HH:MM:SS"
        self.analysis_date = date_format.format(
                                 year = now.year,
                                 month = now.month,
                                 day = now.day,
                                 hour = now.hour,
                                 min = now.minute,
                                 second = now.second)
        # format "YYYYmmdd-HHMMSS"
        self.analysis_timestamp = timestamp_format.format(
                                      year = now.year,
                                      month = now.month,
                                      day = now.day,
                                      hour = now.hour,
                                      min = now.minute,
                                      second = now.second)

        # pipeline version
        #proc = subprocess.Popen(['gcat_workflow_cloud --version 2>&1'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        #self.pipeline_version = (proc.communicate()[0]).split("\n")[0]
        self.pipeline_version = ini.__version__
        
        # path to upload
        [name, ext] = os.path.splitext(os.path.basename(sample_conf_file))
        self.sample_conf_storage_path = "%s/config/%s-%s%s" % (output_dir.rstrip("/"), name, self.analysis_timestamp, ext)
        [name, ext] = os.path.splitext(os.path.basename(param_conf_file))
        self.param_conf_storage_path = "%s/config/%s-%s%s" % (output_dir.rstrip("/"), name, self.analysis_timestamp, ext)
        
    def get_meta_info(self, image):
        #import pwd
        #return "# Docker Image: %s\n# Analysis Date: %s\n# User: %s" % (
        #            image, 
        #            self.analysis_date, 
        #            pwd.getpwuid(os.getuid()).pw_name
        #        )
        return "# Docker Image: %s" % (
                    image
                )

    def get_owner_info(self):
        import pwd
        return pwd.getpwuid(os.getuid()).pw_name

