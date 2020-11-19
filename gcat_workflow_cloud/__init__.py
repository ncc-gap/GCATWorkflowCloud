__version__ = '0.1.1'

import argparse
import gcat_workflow_cloud.run

def main():

    aparser = argparse.ArgumentParser(prog = "gcat_workflow_cloud")
    aparser.add_argument("--version", action = "version", version = "gcat_workflow_cloud-" + __version__ )
    aparser.add_argument('analysis_type', choices=['germline', 'rna'], help = "analysis type")
    aparser.add_argument("sample_conf_file", default = None, type = str, help = "Sample config file")
    aparser.add_argument("output_dir", default = None, type = str,
                         help = "Output directory for Google Cloud Storage or AWS S3 bucket")
    aparser.add_argument("param_conf_file", help = "Parameter config file", type = str)
    aparser.add_argument("--engine", choices = ["awsub", "dsub", "azmon", "ecsub"], default = "ecsub", type = str,
                         help = "Batch job engine")
    aparser.add_argument("--dryrun", help = "For dry run ", action = 'store_true')
    aparser.add_argument("--work_dir", default = "", type = str, help = "The directory path for local work")
    args = aparser.parse_args()

    gcat_workflow_cloud.run.run(args)

