[![Build Status](https://travis-ci.org/ncc-ccat-gap/GCATWorkflowCloud.svg?branch=master)](https://travis-ci.org/ncc-ccat-gap/GCATWorkflowCloud)
![Python](https://img.shields.io/badge/python-3.6%20%7C%203.7-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# GCAT-Workflow Cloud

# Installation

## GCAT-Workflow Cloud

```sh
https://github.com/ncc-ccat-gap/GCATWorkflowCloud.git
cd GCATWorkflowCloud
python setup.py build install
```

## Dependency

See, https://github.com/aokad/ecsub

# Quick Start

```sh
export YOUR_BUCKET=s3://aokad-ana-tokyo

gcat_workflow_cloud germline \
  example_conf/germline_sample.csv \
  ${YOUR_BUCKET}/gcat_workflow_cloud_test \
  example_conf/germline_gcat_ecsub.cfg
```

