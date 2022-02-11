"""
This class implements support for a Slurm Workload Manager cluster.

NOTE that this class (currently) only supports runBatchCommands as method to submit job files.
"""

import os
import re
import subprocess
import sys
import uuid
from subprocess import CalledProcessError, PIPE
from typing import List

from .job_template import job_template
from ..Cluster import Cluster


class SLURM(Cluster):

        def __init__(self):
                super().__init__()
                self.sbatch = 'sbatch'
                self.squeue = 'squeue'

        def getRunningJobIDs(self) -> List[str]:
                return [line.split()[0] for line in self.status()]

        def submit(self, job_file: str) -> str:
                # sbatch example output: b'Submitted batch job 13\n'
                try:
                        p = subprocess.run([self.sbatch, job_file], stdout=PIPE, check=True)
                except CalledProcessError as e:
                        sys.exit(f'Failed to submit job: {e}')
                stdout = p.stdout.decode('utf-8')
                print(stdout, end='')
                return stdout.split()[3]

        def status(self) -> List[str]:
                # squeue example output:
                # JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)\n
                #    90       199 PRosetta     root  R       2:39      1 102-smv4\n'
                p = subprocess.run([self.squeue], stdout=PIPE)
                return [line for line in p.stdout.decode('utf-8').split('\n') if
                        re.search('^\s*[0-9]+', line) is not None]

        def writeJobFile(self, job_file: str, commands: List[str], mem: int) -> None:
                cur_job = open(job_file, 'w')
                job_text = job_template.replace('MEM_PER_CPU', str(mem))
                job_text = job_text.replace('CPUS_PER_TASK', len(commands))
                job_text = job_text.replace('JOB_NAME', 'PRosettaC_' + str(uuid.uuid4()))

                # Prepend scheduler params to job file, if the environment variable SCHEDULER_PARAMS is set.
                loaded_params = False
                if Cluster.SCHEDULER_PARAMS is not None:
                        with open(Cluster.SCHEDULER_PARAMS) as f:
                                try:
                                        sched_params = f.readlines()
                                        job_text = job_text.replace('HEADER', ''.join(sched_params))
                                        loaded_params = True
                                except IOError as e:
                                        print(f'Could not read scheduler params from {Cluster.SCHEDULER_PARAMS}: {e}',
                                              file=sys.stderr)
                if not loaded_params:
                        job_text = job_text.replace('HEADER', '')
                cur_job.write(job_text)
                for command in commands:
                        cur_job.write(command + ' &\n')
                cur_job.write('wait\n')
                cur_job.close()

        def runSingle(self, command):
                raise NotImplementedError('This is not implemented for SGE.')

        def runSingleDepend(self, command, depend_jobs_file):
                raise NotImplementedError('This is not implemented for SGE.')

        def runSingleShell(self, command):
                raise NotImplementedError('This is not implemented for SGE.')

        def runBatchJobs(self, dirlist, command, batch_size=12, mem='2000mb'):
                raise NotImplementedError('This is not implemented for SGE.')
