import os
import re
import socket
import subprocess
import sys
from subprocess import Popen, PIPE
from typing import List

import utils

sys.path.append(utils.SCRIPTS_FOL + '/cluster/PBS')
import chemfarm_job_submission as cjob
import chemfarm_job_batch as cbatch

sys.path.append('.')
from ..Cluster import Cluster

PBS_HOME = os.environ.get("PBS_HOME", "")


class PBS(Cluster):
        def getRunningJobIDs(self) -> List[str]:
                """Return list of IDs of currently running jobs."""
                # not sure if i need to use pbs_Home fpr qstat
                p = Popen([os.path.join(PBS_HOME, 'qstat')], stdout=PIPE, stderr=PIPE)
                jobs_running = [line.split()[0] for line in p.stdout.read().decode("utf-8").split('\n') if re.match("[0-9]+\.\w",line) is not None]
                return jobs_running

        def submit(self, job_file: str) -> str:
                """Submit job_file as job to the cluster and return its assigned ID."""
                p = Popen([os.path.join(PBS_HOME, 'qsub'), job_file], stdout=PIPE, stderr=PIPE)
                s = p.stdout.read()
                job_id = s.strip().decode("utf-8")
                return job_id

        def writeJobFile(self, job_file: str, commands: List[str], mem: int) -> None:
                """Write commands and params to job_file in a scheduler-specific format, so that job_file can be submitted
                        to the cluster."""
                cur_job = open(job_file, 'w')

                cur_job.write(cbatch.sendjob_text[0])

                if Cluster.SCHEDULER_PARAMS != None:
                        SCHEDULER_PARAMS = open(Cluster.SCHEDULER_PARAMS, 'r')
                        for l in SCHEDULER_PARAMS:
                                cur_job.write(l)
                        SCHEDULER_PARAMS.close()

                for line in cbatch.sendjob_text[1:4]:
                        cur_job.write(line)

                cur_job.write(str(len(commands)))
                cur_job.write(cbatch.sendjob_text[4])
                cur_job.write(str(mem) + 'mb')
                for line in cbatch.sendjob_text[5:]:
                        cur_job.write(line)
                for command in commands:
                        cur_job.write(command + ' &\n')
                cur_job.write('wait\n')
                cur_job.close()

        def __init__(self):
                self.typ = "CHEM"

        def runSingle(self, command):
                if(self.typ == "WEXAC"):
                        subprocess.call(["bsub", "-u", "/dev/null", "-R", "rusage[mem=1024]", "-q", "new-all.q", "-o", "out", "-e", "err", command])
                if(self.typ == "CHEM"):
                        job_file = "job_submission.sh"
                        cur_job = open(job_file, 'w')
                        for line in cjob.sendjob_text[:7]:
                                cur_job.write(line)
                        cur_job.write(cjob.sendjob_text[7] + '\'' + command + '\'')
                        for line in cjob.sendjob_text[8:]:
                                cur_job.write(line)
                        cur_job.close()
                        p = Popen([os.path.join(PBS_HOME, 'qsub'), job_file], stdout=PIPE, stderr=PIPE)
                        return [line.split()[0] for line in p.stdout.read().split(b'\n') if b'pbs' in line][0]
        def runSingleDepend(self, command, depend_jobs_file):
                with open(depend_jobs_file, 'r') as f:
                        depend_jobs = [line.split()[0] for line in f]
                if(self.typ == "WEXAC"):
                        subprocess.call(["bsub", "-u", "/dev/null", "-R", "rusage[mem=1024]", "-q", "new-all.q", "-o", "out", "-e", "err", command])
                if(self.typ == "CHEM"):
                        job_file = "job_submission.sh"
                        cur_job = open(job_file, 'w')
                        for line in cjob.sendjob_text[:1]:
                                cur_job.write(line)
                        for line in depend_jobs:
                                cur_job.write('#PBS -W depend=afterany:' + ":".join(depend_jobs))
                        for line in cjob.sendjob_text[1:7]:
                                cur_job.write(line)
                        cur_job.write(cjob.sendjob_text[7] + '\'' + command + '\'')
                        for line in cjob.sendjob_text[8:]:
                                cur_job.write(line)
                        cur_job.close()
                        subprocess.call([os.path.join(PBS_HOME, 'qsub'), job_file])

        def runSingleShell(self, command):
                job_file = "job_submission.sh"
                cur_job = open(job_file, 'w')
                for line in cjob.sendjob_text[:7]:
                        cur_job.write(line)
                cur_job.write(command + '\n')
                for line in cjob.sendjob_text[9:]:
                        cur_job.write(line)
                cur_job.close()
                subprocess.call([os.path.join(PBS_HOME, 'qsub'), job_file])

        def runBatchJobs(self, dirlist, command, batch_size=12, mem=2000):
                curr = os.getcwd()
                f = open(dirlist, 'r')
                lines = [line[:-1] for line in f.readlines()]
                f.close()
                jobs = []
                batch_list = [lines[i:i+batch_size] for i in range(0, len(lines), batch_size)]
                for i, batch in enumerate(batch_list):
                        job_file = "job_batch_" + str(i) + ".sh"
                        cur_job = open(job_file, 'w')
                        for line in cbatch.sendjob_text[:4]:
                                cur_job.write(line)
                        cur_job.write(str(len(batch)))
                        cur_job.write(cbatch.sendjob_text[4])
                        cur_job.write(str(mem) + "mb")
                        for line in cbatch.sendjob_text[5:]:
                                cur_job.write(line)
                        for dirname in batch:
                                cur_job.write('cd ' + dirname + '\n')
                                cur_job.write(command + ' &\n')
                                cur_job.write('cd ' + curr + '\n')
                        cur_job.write('wait\n')
                        cur_job.close()
                        p = Popen(["qsub", job_file], stdout=PIPE, stderr=PIPE)
                        jobs.append([line.split()[0] for line in p.stdout.read().split(b'\n') if b'pbs' in line][0])
                return jobs

# Inner functions
        def getServer(self):
                server_name = socket.gethostname()
                if 'chemfarm' in server_name:
                        self.typ = "CHEM"
                elif 'wexac' in server_name:
                        self.typ = "WEXAC"
                else:
                        print("This is not a cluster supporting server")
                        sys.exit()
