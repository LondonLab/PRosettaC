import utils
import os,sys
sys.path.append(utils.SCRIPTS_FOL + 'PBS/')
import cluster_pbs

def main(name, argv):
        if len(argv) != 1:
                print_usage(name)
                return

        cluster = cluster_pbs.cluster_pbs()
        cluster.runSingle("python " + utils.SCRIPTS_FOL + 'main.py ' + argv[0])

def print_usage(name):
        print("Usage : " + name + " <param file>")

if __name__ == "__main__":
    main(sys.argv[0], sys.argv[1:])
