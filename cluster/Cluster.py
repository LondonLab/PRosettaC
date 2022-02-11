import os
import sys
import time
from typing import List


class Cluster:
        """
        Interface between PRosettaC's main algorithm and job scheduler.

        This class defines all necessary methods do communicate with the job scheduler and is used in the main scripts. The
        methods handling scheduler-specific code have to be implemented by a dedicated subclass. In the main scripts,
        cluster.getCluster gives you an object of the correct subclass to talk to the scheduler.
        """

        # Prepend content of this file to each job file, if set.
        SCHEDULER_PARAMS = os.environ.get('SCHEDULER_PARAMS', None)

        def getRunningJobIDs(self) -> List[str]:
                """Return list of IDs of currently running jobs."""
                raise NotImplementedError()

        def submit(self, job_file: str) -> str:
                """Submit job_file as job to the cluster and return its assigned ID."""
                raise NotImplementedError()

        def writeJobFile(self, job_file: str, commands: List[str], mem: int) -> None:
                """Write commands and params to job_file in a scheduler-specific format, so that job_file can be submitted
                to the cluster."""
                raise NotImplementedError()

        def runSingle(self, command):
                raise NotImplementedError()

        def runSingleDepend(self, command, depend_jobs_file):
                raise NotImplementedError()

        def runSingleShell(self, command):
                raise NotImplementedError()

        def runBatchJobs(self, dirlist, command, batch_size=12, mem='2000mb'):
                raise NotImplementedError()

        def wait(self, job_ids: List[str], timeout=-1) -> None:
                """
                Wait until all jobs in job_ids have finished or the optional timeout has been reached.
                :param job_ids: list of job IDs to wait for
                :type job_ids: List[str]
                :param timeout: number of seconds (in thousands) before timeout (default none)
                :type timeout: int
                """
                i = 0
                while self.jobsRunning(job_ids):
                        print("Not done yet")
                        if i == timeout:
                                break
                        time.sleep(1)
                        i += 1

        def jobsRunning(self, job_ids: List[str]) -> bool:
                """Return true, if any job in job_ids is currently running, false otherwise."""
                jobs_running = self.getRunningJobIDs()
                for id in job_ids:
                        if id in jobs_running:
                                return True
                else:
                        return False

        def runBatchCommands(self, commands: List[str], mem, batch_size=12) -> List[str]:
                """
                Submit commands as jobs to the cluster in batches of batch_size. Return the assigned job IDs.
                :param commands: Each element is a command placed on its own line in a job file
                :type commands: List[str]
                :param batch_size: numbers of commands per job file (default 12)
                :type batch_size: int
                :param mem: amount of memory assigned to each job (in MB)
                :type params: int
                :return: list of IDs of submitted jobs
                :rtype: List[str]
                """
                jobs = []
                batch_list = [commands[i:i + batch_size] for i in range(0, len(commands), batch_size)]
                for i, batch in enumerate(batch_list):
                        job_file = "job_batch_" + str(i) + ".sh"
                        self.writeJobFile(job_file, commands=batch, mem=mem)
                        jobs.append(self.submit(job_file))
                return jobs

        def runCommands(self, commands):
                jobs = []
                for command in commands:
                        jobs.append(self.runSingle(command))
                return jobs

        def runDirSingle(self, dirname, command):
                pyutils.create_folder(dirname[:-1])
                curr = os.getcwd()
                os.chdir(dirname[:-1])
                job = self.runSingle(command)
                os.chdir(curr)
                return job

        def runDirCommands(self, dirlist, commands):
                self.checkList(dirlist, commands)
                f = open(dirlist, 'r')
                lines = f.readlines()
                jobs = []
                for dirname, command in zip(lines, commands):
                        jobs.append(self.runDirSingle(dirname, command))
                f.close()
                return jobs

        def runJobs(self, dirlist, command):
                f = open(dirlist, 'r')
                lines = f.readlines()
                commands = [command] * len(lines)
                f.close()
                return self.runDirCommands(dirlist, commands)

        def runJobsArgs(self, dirlist, command, arglist):
                self.checkList(dirlist, arglist)
                commands = []
                for arg in arglist:
                        commands += [command + ' ' + str(arg)]
                self.runDirCommands(dirlist, commands)

        def runCommandsArgs(self, command, arglist):
                commands = []
                for arg in arglist:
                        commands.append(command + ' ' + str(arg))
                return self.runCommands(commands)

        def runJobsName(self, dirlist, command):
                with open(dirlist) as f:
                        lines = f.read().splitlines()
                self.runJobsArgs(dirlist, command, lines)

        def checkList(self, dirlist, l):
                f = open(dirlist, 'r')
                if not len(f.readlines()) == len(l):
                        print('length of list does not match dirlist length')
                        sys.exit()
                f.close()
