#! /usr/bin/env python

import os, multiprocessing, configparser
import gcat_workflow_cloud.batch_engine as be
import gcat_workflow_cloud.run_conf as rc
import gcat_workflow_cloud.storage as st

class JobError(Exception):
    def error_text(name, code):
        return "[%s] Failure job execution with exit_code (%s)" % (name, str(code))

def run(args):

    if args.analysis_type == "germline":
        import gcat_workflow_cloud.sample_conf_germline as sc
    elif args.analysis_type == "somatic":
        import gcat_workflow_cloud.sample_conf_somatic as sc
    elif args.analysis_type == "rna":
        import gcat_workflow_cloud.sample_conf_rna as sc

    args.output_dir = args.output_dir.rstrip("/")
    
    sample_conf = sc.Sample_conf(args.sample_conf_file, exist_check = False)
    
    param_conf = configparser.ConfigParser()
    param_conf.read(args.param_conf_file)

    run_conf = rc.Run_conf(
        sample_conf_file = args.sample_conf_file, 
        param_conf_file = args.param_conf_file,
        analysis_type = args.analysis_type,
        output_dir = args.output_dir,
        use_bam = args.use_bam
    )
                        
    # temporary procedure
    if args.work_dir == "":
        tmp_dir = "%s/tmp/%s" % (os.getcwd(), run_conf.analysis_timestamp)
    else:
        tmp_dir = "%s/%s" % (args.work_dir, run_conf.analysis_timestamp)

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
        print ("Creating temporary directory: " +  tmp_dir)

    # preparing batch job engine
    if args.engine == "dsub":
        batch_engine = be.Dsub_factory()
    elif args.engine == "azmon":
        batch_engine = be.Azmon_factory()
    elif args.engine == "ecsub":
        batch_engine = be.Ecsub_factory()
        batch_engine.s3_wdir = run_conf.output_dir + "/ecsub"
        batch_engine.wdir = tmp_dir + "/ecsub"
    else:
        batch_engine = be.Awsub_factory()

    batch_engine.dryrun = args.dryrun
    batch_engine.general_param = param_conf.get("general", "instance_option")
    
    # upload config files
    storage = st.Storage(dryrun = args.dryrun)
    storage.upload(args.sample_conf_file, run_conf.sample_conf_storage_path, create_bucket = True)
    storage.upload(args.param_conf_file, run_conf.param_conf_storage_path, create_bucket = True)

    ##########
    # germline
    if args.analysis_type == "germline":
        
        # download metadata files
        sample_conf.readgroup_local = {}
        readgroup_dir = tmp_dir + "/readgroup"
        os.makedirs(readgroup_dir, exist_ok=True)
        for sample in sample_conf.readgroup:
            local_path = "%s/%s.txt" % (readgroup_dir, sample)
            storage.download(local_path, sample_conf.readgroup[sample])
            sample_conf.readgroup_local[sample] = local_path
        
        import gcat_workflow_cloud.tasks.bamtofastq as bamtofastq
        import gcat_workflow_cloud.tasks.gatk_fq2cram as fq2cram
        import gcat_workflow_cloud.tasks.gatk_haploypecaller as haploypecaller
        import gcat_workflow_cloud.tasks.gatk_collectwgsmetrics as collectwgsmetrics
        import gcat_workflow_cloud.tasks.gatk_collectmultiplemetrics as collectmultiplemetrics
        import gcat_workflow_cloud.tasks.gridss as gridss
        import gcat_workflow_cloud.tasks.manta as manta
        import gcat_workflow_cloud.tasks.melt as melt
        import gcat_workflow_cloud.tasks.fastqc as fastqc

        bamtofastq_task = bamtofastq.Task(tmp_dir, sample_conf, param_conf, run_conf)
        fq2cram_task = fq2cram.Task(tmp_dir, sample_conf, param_conf, run_conf)
        haploypecaller_task_list = []
        for key in sorted(haploypecaller.TAGS.keys()):
            haploypecaller_task_list.append(haploypecaller.Task(tmp_dir, sample_conf, param_conf, run_conf, key))
        collectwgsmetrics_task = collectwgsmetrics.Task(tmp_dir, sample_conf, param_conf, run_conf)
        collectmultiplemetrics_task = collectmultiplemetrics.Task(tmp_dir, sample_conf, param_conf, run_conf)
        gridss_task = gridss.Task(tmp_dir, sample_conf, param_conf, run_conf)
        manta_task = manta.Task(tmp_dir, sample_conf, param_conf, run_conf)
        melt_task = melt.Task(tmp_dir, sample_conf, param_conf, run_conf)
        fastqc_task = fastqc.Task(tmp_dir, sample_conf, param_conf, run_conf)

        # bamtofastq
        p_bamtofastq = multiprocessing.Process(target = batch_engine.execute, args = (bamtofastq_task,))
        try:
            p_bamtofastq.start()
            p_bamtofastq.join()
        except KeyboardInterrupt:
            pass
        
        if p_bamtofastq.exitcode != 0:
            raise JobError(JobError.error_text(bamtofastq_task.TASK_NAME, p_bamtofastq.exitcode))

        # fq2cram
        p_fq2cram = multiprocessing.Process(target = batch_engine.execute, args = (fq2cram_task,))
        p_fastqc = multiprocessing.Process(target = batch_engine.execute, args = (fastqc_task,))

        try:
            p_fq2cram.start()
            p_fastqc.start()
            p_fq2cram.join()
        except KeyboardInterrupt:
            pass
        
        if p_fq2cram.exitcode != 0:
            raise JobError(JobError.error_text(fq2cram_task.TASK_NAME, p_fq2cram.exitcode))

        # analysis
        p_haploypecaller = multiprocessing.Process(target = batch_engine.seq_execute, args = (haploypecaller_task_list, ))
        p_collectwgsmetrics = multiprocessing.Process(target = batch_engine.execute, args = (collectwgsmetrics_task,))
        p_collectmultiplemetrics = multiprocessing.Process(target = batch_engine.execute, args = (collectmultiplemetrics_task,))
        p_gridss = multiprocessing.Process(target = batch_engine.execute, args = (gridss_task,))
        p_manta = multiprocessing.Process(target = batch_engine.execute, args = (manta_task,))
        p_melt = multiprocessing.Process(target = batch_engine.execute, args = (melt_task,))

        try:
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
            p_fastqc.join()

        except KeyboardInterrupt:
            pass

        if p_haploypecaller.exitcode != 0:
            raise JobError(JobError.error_text(haploypecaller_task_list[p_haploypecaller.exitcode-1].TASK_NAME, p_haploypecaller.exitcode))
        
        if p_collectwgsmetrics.exitcode != 0:
            raise JobError(JobError.error_text(collectwgsmetrics_task.TASK_NAME, p_collectwgsmetrics.exitcode))
        
        if p_collectmultiplemetrics.exitcode != 0:
            raise JobError(JobError.error_text(collectmultiplemetrics_task.TASK_NAME, p_collectmultiplemetrics.exitcode))
        
        if p_gridss.exitcode != 0:
            raise JobError(JobError.error_text(gridss_task.TASK_NAME, p_gridss.exitcode))
        
        if p_manta.exitcode != 0:
            raise JobError(JobError.error_text(manta_task.TASK_NAME, p_manta.exitcode))
        
        if p_melt.exitcode != 0:
            raise JobError(JobError.error_text(melt_task.TASK_NAME, p_melt.exitcode))
        
        if p_fastqc.exitcode != 0:
            raise JobError(JobError.error_text(fastqc_task.TASK_NAME, p_fastqc.exitcode))

    ##########
    # somatic
    if args.analysis_type == "somatic":
        
        # download metadata files
        sample_conf.readgroup_local = {}
        readgroup_dir = tmp_dir + "/readgroup"
        os.makedirs(readgroup_dir, exist_ok=True)
        for sample in sample_conf.readgroup:
            local_path = "%s/%s.txt" % (readgroup_dir, sample)
            storage.download(local_path, sample_conf.readgroup[sample])
            sample_conf.readgroup_local[sample] = local_path
        
        import gcat_workflow_cloud.tasks.bamtofastq as bamtofastq
        import gcat_workflow_cloud.tasks.gatk_fq2cram as fq2cram
        import gcat_workflow_cloud.tasks.gatk_mutectcaller as mutectcaller
        import gcat_workflow_cloud.tasks.gatk_collectwgsmetrics as collectwgsmetrics
        import gcat_workflow_cloud.tasks.gatk_collectmultiplemetrics as collectmultiplemetrics
        import gcat_workflow_cloud.tasks.gridss as gridss
        import gcat_workflow_cloud.tasks.manta as manta
        import gcat_workflow_cloud.tasks.genomonsv as genomonsv
        import gcat_workflow_cloud.tasks.genomon_mutation_call as genomon_mutation_call
        import gcat_workflow_cloud.tasks.fastqc as fastqc

        bamtofastq_task = bamtofastq.Task(tmp_dir, sample_conf, param_conf, run_conf)
        fq2cram_task = fq2cram.Task(tmp_dir, sample_conf, param_conf, run_conf)
        mutectcaller_task = mutectcaller.Task(tmp_dir, sample_conf, param_conf, run_conf)
        collectwgsmetrics_task = collectwgsmetrics.Task(tmp_dir, sample_conf, param_conf, run_conf)
        collectmultiplemetrics_task = collectmultiplemetrics.Task(tmp_dir, sample_conf, param_conf, run_conf)
        gridss_task = gridss.Task(tmp_dir, sample_conf, param_conf, run_conf)
        manta_task = manta.Task(tmp_dir, sample_conf, param_conf, run_conf)
        genomonsv_task = genomonsv.Task(tmp_dir, sample_conf, param_conf, run_conf)
        genomon_mutation_call_task = genomon_mutation_call.Task(tmp_dir, sample_conf, param_conf, run_conf)
        fastqc_task = fastqc.Task(tmp_dir, sample_conf, param_conf, run_conf)

        # bamtofastq
        p_bamtofastq = multiprocessing.Process(target = batch_engine.execute, args = (bamtofastq_task,))
        try:
            p_bamtofastq.start()
            p_bamtofastq.join()
        except KeyboardInterrupt:
            pass
        
        if p_bamtofastq.exitcode != 0:
            raise JobError(JobError.error_text(bamtofastq_task.TASK_NAME, p_bamtofastq.exitcode))

        # fq2cram
        p_fq2cram = multiprocessing.Process(target = batch_engine.execute, args = (fq2cram_task,))
        p_fastqc = multiprocessing.Process(target = batch_engine.execute, args = (fastqc_task,))

        try:
            p_fq2cram.start()
            p_fastqc.start()
            p_fq2cram.join()
        except KeyboardInterrupt:
            pass
        
        if p_fq2cram.exitcode != 0:
            raise JobError(JobError.error_text(fq2cram_task.TASK_NAME, p_fq2cram.exitcode))

        # analysis
        p_mutectcaller = multiprocessing.Process(target = batch_engine.execute, args = (mutectcaller_task, ))
        p_collectwgsmetrics = multiprocessing.Process(target = batch_engine.execute, args = (collectwgsmetrics_task,))
        p_collectmultiplemetrics = multiprocessing.Process(target = batch_engine.execute, args = (collectmultiplemetrics_task,))
        p_gridss = multiprocessing.Process(target = batch_engine.execute, args = (gridss_task,))
        p_manta = multiprocessing.Process(target = batch_engine.execute, args = (manta_task,))
        p_genomonsv = multiprocessing.Process(target = batch_engine.execute, args = (genomonsv_task,))
        p_genomon_mutation_call = multiprocessing.Process(target = batch_engine.execute, args = (genomon_mutation_call_task,))

        try:
            p_mutectcaller.start()
            p_collectwgsmetrics.start()
            p_collectmultiplemetrics.start()
            p_gridss.start()
            p_manta.start()
            p_genomonsv.start()
            p_genomon_mutation_call.start()

            p_mutectcaller.join()
            p_collectwgsmetrics.join()
            p_collectmultiplemetrics.join()
            p_gridss.join()
            p_manta.join()
            p_genomonsv.join()
            p_genomon_mutation_call.join()
            p_fastqc.join()

        except KeyboardInterrupt:
            pass

        if p_mutectcaller.exitcode != 0:
            raise JobError(JobError.error_text(mutectcaller_task.TASK_NAME, p_mutectcaller.exitcode))
        
        if p_collectwgsmetrics.exitcode != 0:
            raise JobError(JobError.error_text(collectwgsmetrics_task.TASK_NAME, p_collectwgsmetrics.exitcode))
        
        if p_collectmultiplemetrics.exitcode != 0:
            raise JobError(JobError.error_text(collectmultiplemetrics_task.TASK_NAME, p_collectmultiplemetrics.exitcode))
        
        if p_gridss.exitcode != 0:
            raise JobError(JobError.error_text(gridss_task.TASK_NAME, p_gridss.exitcode))
        
        if p_manta.exitcode != 0:
            raise JobError(JobError.error_text(manta_task.TASK_NAME, p_manta.exitcode))
        
        if p_genomonsv.exitcode != 0:
            raise JobError(JobError.error_text(genomonsv_task.TASK_NAME, p_genomonsv.exitcode))
        
        if p_genomon_mutation_call.exitcode != 0:
            raise JobError(JobError.error_text(genomon_mutation_call_task.TASK_NAME, p_genomon_mutation_call.exitcode))
        
        if p_fastqc.exitcode != 0:
            raise JobError(JobError.error_text(fastqc_task.TASK_NAME, p_fastqc.exitcode))

    ##########
    # rna
    elif args.analysis_type == "rna":

        import gcat_workflow_cloud.tasks.star_align as star_align
        import gcat_workflow_cloud.tasks.expression as expression
        import gcat_workflow_cloud.tasks.intron_retention as intron_retention

        star_align_task = star_align.Task(tmp_dir, sample_conf, param_conf, run_conf)
        expression_task = expression.Task(tmp_dir, sample_conf, param_conf, run_conf)
        intron_retention_task = intron_retention.Task(tmp_dir, sample_conf, param_conf, run_conf)

        p_star_align = multiprocessing.Process(target = batch_engine.execute, args = (star_align_task,))

        try:
            p_star_align.start()
            p_star_align.join()
        except KeyboardInterrupt:
            pass
        
        if p_star_align.exitcode != 0:
            raise JobError(JobError.error_text(star_align_task.TASK_NAME, p_star_align.exitcode))

        p_expression = multiprocessing.Process(target = batch_engine.execute, args = (expression_task,))
        p_intron_retention = multiprocessing.Process(target = batch_engine.execute, args = (intron_retention_task,))

        try:
            p_expression.start()
            p_intron_retention.start()

            p_expression.join()
            p_intron_retention.join()

        except KeyboardInterrupt:
            pass
        
        if p_expression.exitcode != 0:
            raise JobError(JobError.error_text(expression_task.TASK_NAME, p_expression.exitcode))
        
        if p_intron_retention.exitcode != 0:
            raise JobError(JobError.error_text(intron_retention_task.TASK_NAME, p_intron_retention.exitcode))
        
    if args.engine == "ecsub":
        if args.skip_price == True:
            print ("Skip print price.")
        else:
            metrics_dir = tmp_dir + "/metrics"
            os.makedirs(metrics_dir, exist_ok=True)
            
            metrics_file = metrics_dir + "/summary.json"
            files = batch_engine.print_metrics(metrics_file)
            storage.upload(metrics_file, run_conf.output_dir + "/metrics/metrics.json")
            for f in files:
                storage.upload(f, "%s/%s" % (batch_engine.s3_wdir, "/".join(f.split("/")[-3:])))

            cost_file = metrics_dir + "/cost.csv"
            batch_engine.print_summary(cost_file)
            storage.upload(cost_file, run_conf.output_dir + "/metrics/cost.csv")
