import utils
import os,sys
sys.path.append(utils.SCRIPTS_FOL)
import cluster as cl

def main(name, argv):
        if len(argv) != 1:
                print_usage(name)
                return

        params = utils.read_params(argv[0])
        cluster = cl.getCluster(params['ClusterName'])
        cluster.runSingle("python " + utils.SCRIPTS_FOL + 'auto.py ' + argv[0])

def print_usage(name):
        print("Usage : " + name + " <param file>")

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
