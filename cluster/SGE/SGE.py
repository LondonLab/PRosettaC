import sys
from typing import List
import re
import subprocess
from subprocess import CalledProcessError

sys.path.append('.')
from job_template import job_template
from ..Cluster import Cluster


class SGE(Cluster):
    def getRunningJobIDs(self) -> List[str]:
        pass

    def submit(self, job_file: str) -> str:
        pass

    def status(self) -> List[str]:
        pass

    def writeJobFile(self, job_file: str, commands: List[str], params: dict) -> None:
        pass
