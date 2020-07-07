#! /usr/bin/env python

import sys, os, multiprocessing
   
import gcat_workflow_cloud.batch_engine as be
import gcat_workflow_cloud.sample_conf as sc
import gcat_workflow_cloud.run_conf as rc
import gcat_workflow_cloud.storage as st

if sys.version_info.major == 2:
    import ConfigParser as cp
else:
    import configparser as cp

class JobError(Exception):
    template = "[%s] Failure job execution with exit_code (%d)"

def run(args):
    args.output_dir = args.output_dir.rstrip("/")
    
    sample_conf = sc.Sample_conf(args.sample_conf_file, exist_check = False)
    
    param_conf = cp.ConfigParser()
    param_conf.read(args.param_conf_file)

    run_conf = rc.Run_conf(sample_conf_file = args.sample_conf_file, 
                        param_conf_file = args.param_conf_file,
                        analysis_type = args.analysis_type,
                        output_dir = args.output_dir)
                        
    # tmp_dir = tempfile.mkdtemp()
    # temporary procedure
    tmp_dir = os.getcwd() + "/tmp"
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
        print ("Creating temporary directory: " +  tmp_dir)

    log_dir = os.getcwd() + "/log"
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
        print ("Creating log directory: " +  log_dir)
        
    # preparing batch job engine
    if args.engine == "dsub":
        factory = be.Dsub_factory()
    elif args.engine == "azmon":
        factory = be.Azmon_factory()
    elif args.engine == "ecsub":
        factory = be.Ecsub_factory()
        factory.s3_wdir = args.output_dir + "/ecsub"
        factory.wdir = tmp_dir + "/ecsub"
    else:
        factory = be.Awsub_factory()

    factory.dryrun = args.dryrun
    batch_engine = be.Batch_engine(factory, param_conf.get("general", "instance_option"))
    
    # upload config files
    storage = st.Storage(dryrun = args.dryrun)
    storage.upload(args.sample_conf_file, run_conf.sample_conf_storage_path, create_bucket = True)
    storage.upload(args.param_conf_file, run_conf.param_conf_storage_path, create_bucket = True)
    
    ##########
    # germline
    if args.analysis_type == "germline":

        import gcat_workflow_cloud.tasks.gatk_fq2cram as fq2cram
        import gcat_workflow_cloud.tasks.gatk_haploypecaller as haploypecaller
        import gcat_workflow_cloud.tasks.gatk_collectwgsmetrics as collectwgsmetrics
        import gcat_workflow_cloud.tasks.gatk_collectmultiplemetrics as collectmultiplemetrics
        import gcat_workflow_cloud.tasks.gridss as gridss
        import gcat_workflow_cloud.tasks.manta as manta
        import gcat_workflow_cloud.tasks.melt as melt

        fq2cram_task = fq2cram.Task(args.output_dir, tmp_dir, sample_conf, param_conf, run_conf)
        haploypecaller_task = haploypecaller.Task(args.output_dir, tmp_dir, sample_conf, param_conf, run_conf)
        collectwgsmetrics_task = collectwgsmetrics.Task(args.output_dir, tmp_dir, sample_conf, param_conf, run_conf)
        collectmultiplemetrics_task = collectmultiplemetrics.Task(args.output_dir, tmp_dir, sample_conf, param_conf, run_conf)
        gridss_task = gridss.Task(args.output_dir, tmp_dir, sample_conf, param_conf, run_conf)
        manta_task = manta.Task(args.output_dir, tmp_dir, sample_conf, param_conf, run_conf)
        melt_task = melt.Task(args.output_dir, tmp_dir, sample_conf, param_conf, run_conf)

        p_fq2cram = multiprocessing.Process(target = batch_engine.execute, args = (fq2cram_task,))
        p_fq2cram.start()
        p_fq2cram.join()
        if p_fq2cram.exitcode != 0:
            raise JobError(JobError.template % (fq2cram_task.TASK_NAME, p_fq2cram.exitcode))

        p_haploypecaller = multiprocessing.Process(target = batch_engine.execute, args = (haploypecaller_task,))
        p_collectwgsmetrics = multiprocessing.Process(target = batch_engine.execute, args = (collectwgsmetrics_task,))
        p_collectmultiplemetrics = multiprocessing.Process(target = batch_engine.execute, args = (collectmultiplemetrics_task,))
        p_gridss = multiprocessing.Process(target = batch_engine.execute, args = (gridss_task,))
        p_manta = multiprocessing.Process(target = batch_engine.execute, args = (manta_task,))
        p_melt = multiprocessing.Process(target = batch_engine.execute, args = (melt_task,))

        p_haploypecaller.start()
        p_collectwgsmetrics.start()
        p_collectmultiplemetrics.start()
        p_gridss.start()
        p_manta.start()
        p_melt.start()

        p_haploypecaller.join()
        p_collectwgsmetrics.join()
        p_collectmultiplemetrics.join()
        p_gridss.join()
        p_manta.join()
        p_melt.join()
    
        if p_fq2cram.exitcode != 0:
            raise JobError(JobError.template % (haploypecaller_task.TASK_NAME, p_haploypecaller.exitcode))
        
        if p_fq2cram.exitcode != 0:
            raise JobError(JobError.template % (collectwgsmetrics_task.TASK_NAME, p_collectwgsmetrics.exitcode))
        
        if p_fq2cram.exitcode != 0:
            raise JobError(JobError.template % (collectmultiplemetrics_task.TASK_NAME, p_collectmultiplemetrics.exitcode))
        
        if p_fq2cram.exitcode != 0:
            raise JobError(JobError.template % (gridss_task.TASK_NAME, p_gridss.exitcode))
        
        if p_fq2cram.exitcode != 0:
            raise JobError(JobError.template % (manta_task.TASK_NAME, p_manta.exitcode))
        
        if p_fq2cram.exitcode != 0:
            raise JobError(JobError.template % (melt_task.TASK_NAME, p_melt.exitcode))
        
    if args.engine == "ecsub":
        factory.print_summary(run_conf, log_dir)
