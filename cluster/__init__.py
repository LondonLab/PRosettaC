import sys

from .Cluster import Cluster
from .SGE.SGE import SGE


def getCluster(name: str) -> Cluster:
        """Get a handle to the cluster given by name."""
        if name == 'SGE':
                return SGE()
        else:
                sys.exit(f'Queue "{name}" not supported.')
