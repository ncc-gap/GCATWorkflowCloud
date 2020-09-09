#!/usr/bin/env python

from setuptools import setup, find_packages
from gcat_workflow_cloud import __version__

setup(
    name = 'gcat_workflow_cloud',
    version = __version__,
    description = 'Python tools for executing GCAT-Workflow in cloud environments.',
    url = 'https://github.com/ncc-ccat-gap/GCATWorkflowCloud.git',
    author = 'Yuichi Shiraishi, Ai Okada, Kenichi Chiba',
    author_email = 'genomon.devel@gmail.com',
    license = 'GPLv3',
    
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Genome Analysis :: RNA-seq',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],

    packages = find_packages(exclude = ['docker']),
    package_data={'gcat_workflow_cloud': ['script/*']},
    # install_requires = ['dsub'],

    entry_points = {'console_scripts': ['gcat_workflow_cloud = gcat_workflow_cloud:main']},
    test_suite = 'tests'

)

