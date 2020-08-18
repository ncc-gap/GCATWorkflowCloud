#! /usr/bin/env python

import subprocess
import signal

class Abstract_factory(object):
    def __init__(self):
        self.general_param = ""
        self.dryrun = False
        
    #@abc.abstractmethod
    def generate_commands(self, task):
        pass
    
    def _execute(self, task):
        proc = subprocess.Popen(self.generate_commands(task))
        
        try:
            while proc.poll() is None:
                pass
            
        except KeyboardInterrupt:
            proc.send_signal(signal.SIGINT)
            while proc.poll() is None:
                pass
            return 1
        
        if proc.returncode == 0:
            return 0
        return 1
        
    def execute(self, task):
        exit(self._execute(task))

    def seq_execute(self, tasks):
        for task in tasks:
            if self._execute(task) != 0:
                exit(1)
        exit(0)
        
    def base_commands(self, commands):
        if self.dryrun:
            return ["echo", ">>"] + commands
        return commands

class Dsub_factory(Abstract_factory):

    def generate_commands(self, task):

        commands = ["dsub"] + self.general_param.split() + task.resource_param.split() + \
                     ["--logging", task.log_dir, "--script", task.script_file, \
                      "--image", task.image, "--tasks", task.task_file, "--wait"]
        return self.base_commands(commands)

class Azmon_factory(Abstract_factory):

    def generate_commands(self, task):

        commands = ["azurebatchmon"] + self.general_param.split() + task.resource_param.split() + \
                     ["--script", task.script_file, "--image", task.image, "--tasks", task.task_file]

        return self.base_commands(commands)

class Awsub_factory(Abstract_factory):

    def generate_commands(self, task):

        commands = ["awsub"] + self.general_param.split() + task.resource_param.split() + \
                     ["--script", task.script_file, "--image", task.image, "--tasks", task.task_file]

        return self.base_commands(commands)

class Ecsub_factory(Abstract_factory):
    
    s3_wdir = ""
    wdir = "/tmp/ecsub"

    def generate_commands(self, task):
        
        commands = ["ecsub", "submit"] + self.general_param.split() + task.resource_param.split() + \
                     ["--script", task.script_file, "--image", task.image, "--tasks", task.task_file] + \
                     ["--aws-s3-bucket", self.s3_wdir, "--wdir", self.wdir]
        import re
        if re.match(r"[0-9]+\.dkr\.ecr\..+\.amazonaws\.com", task.image) != None:
            commands.append("--use-amazon-ecr")

        return self.base_commands(commands)
    
    def print_summary(self, log_file):
        def accessor(obj, key):
            if key in obj:
                return obj[key]
            t_key_list = [s for s in obj.keys() if key.lower().replace("_", "").replace("-", "").replace(".", "") == s.lower().replace("_", "").replace("-", "").replace(".", "")]
            if len(t_key_list) == 0:
                return None
            return obj[t_key_list[0]]
        
        import glob
        import json
        
        dirs = glob.glob(self.wdir + "/*")
        
        sim_cost = 0.0
        od_cost = 0.0
        jobs = []
        header = ""
        for d in sorted(dirs):
            files = glob.glob("%s/log/summary*.log" % (d))
            for f in sorted(files):
                summary = json.load(open(f))
                for j in accessor(summary, "Jobs"):
                    job = {
                        "ClusterName": accessor(summary, "ClusterName"),
                        "Ec2InstanceType": accessor(j, "Ec2InstanceType"),
                        "Ec2InstanceDiskSize": accessor(summary, "Ec2InstanceDiskSize"),
                        "Spot": accessor(j, "Spot"),
                        "WorkHours": accessor(j, "WorkHours")
                    }
                    job["OdPrice"] = accessor(j, "OdPrice")
                    if accessor(j, "Spot"):
                        job["SpotPrice"] = accessor(j, "SpotPrice")
                        job["WorkHours*Price"] = accessor(j, "WorkHours") * accessor(j, "SpotPrice")
                    else:
                        job["WorkHours*Price"] = accessor(j, "WorkHours") * accessor(j, "OdPrice")
                    
                    if len(header) == 0:
                        header = ",".join(sorted(job.keys()))
                        
                    jobs.append(job)
                    sim_cost += job["WorkHours*Price"]
                    od_cost += accessor(j, "WorkHours") * accessor(j, "OdPrice")
        
        print ("Total cost of this instance usage is %.3f USD. See detail, %s" % (sim_cost, log_file))
        
        t = "TotalInstanceCost,%.3f\n" % (sim_cost)
        t += "OndemandCost,%.3f\n" % (od_cost)
        t += "%s\n" % (header)
        for j in jobs:
            for k in sorted(j.keys()):
                if type(j[k]) == type(0.0):
                    t += "%.3f," % (j[k])
                else:
                    t += "%s," % (j[k])
            t += "\n"
        open(log_file, "w").write(t)

    def print_metrics(self, log_file):
        
        import glob
        import json
        
        files = glob.glob(self.wdir + "/*/metrics/*.txt")
        
        jobs = {}
        header = ""
        for f in sorted(files):
            header = True
            data = []
            
            for row in open(f).readlines():
                if header:
                    header = False
                    continue
                if row.rstrip() == "":
                    continue
                data.append(row.split("\t")[2])

            (idx, metrics_name) = f.split("/")[-1].replace(".txt", "").split("-")
            job_name = f.split("/")[-3] + "." + idx
            if not job_name in jobs:
                jobs[job_name] = {}
            jobs[job_name][metrics_name] = max(data)
        
        if len(jobs) > 0:
            metrics_keys = sorted(list(jobs[list(jobs.keys())[0]].keys()))
            print("MaxMetrics (%)")
            print("\t".join(metrics_keys) + "\tJobName")
            for j in jobs:
                text = []
                for m in metrics_keys:
                    text.append(jobs[j][m])
                text.append(j)
                print("\t".join(text))

        json.dump(jobs, open(log_file, "w"), indent=4, sort_keys=True, separators=(',', ': '))
        return files

if __name__ == "__main__":
    pass
