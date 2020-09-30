"""
This class implements support for a Sun Grid Engine cluster.

NOTE that this class (currently) only supports runBatchCommands as method to submit job files.
"""

import os
import re
import subprocess
import sys
from subprocess import CalledProcessError, PIPE
from typing import List

from .job_template import job_template
from ..Cluster import Cluster


class SGE(Cluster):

        def __init__(self):
                super().__init__()
                SGE_HOME = os.environ.get('SGE_HOME', '')
                self.qsub = os.path.join(SGE_HOME, 'qsub')
                self.qstat = os.path.join(SGE_HOME, 'qstat')

        def getRunningJobIDs(self) -> List[str]:
                return [line.split()[0] for line in self.status()]

        def submit(self, job_file: str) -> str:
                # qsub example output: b'Your job 25658 ("job_batch_0.sh") has been submitted\n'
                try:
                        p = subprocess.run([self.qsub, job_file], stdout=PIPE, check=True)
                except CalledProcessError as e:
                        sys.exit(f'Failed to submit job: {e}')
                stdout = p.stdout.decode('utf-8')
                print(stdout, end='')
                return stdout.split()[2]

        def status(self) -> List[str]:
                # qstat example output:
                # b'job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID\n
                # -----------------------------------------------------------------------------------------------------------------\n
                # 26896 2.00000 queuetest. praktikum01  r     08/31/2020 19:28:38 all.q@gaia                         1'
                p = subprocess.run([self.qstat], stdout=PIPE)
                return [line for line in p.stdout.decode('utf-8').split('\n') if
                        re.search('^\s*[0-9]+', line) is not None]

        def writeJobFile(self, job_file: str, commands: List[str], mem: int) -> None:
                cur_job = open(job_file, 'w')
                job_text = job_template.replace('MEMORY', str(mem) + 'M')

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
