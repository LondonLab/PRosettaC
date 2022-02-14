#This is a script for post-processing sub-clustering of individual clusters
#it can help choosing a sinlge model on which to work on
#run from a folder of a specific cluster, it generates a Results/ folder
#with the sub-clusteing.
#It receives only one argument, the RMSD for the sub-clustering
#the default PRosettaC clustering is 4A RMSD, so 4 or above is useless.
#Also, it is only useful for big clusters.
#You can do several rounds, depends on the size of your cluster.
#For example, 2A on the top cluster, then 1A on the top sub-cluster,
#and if still big enough, even a third 0.5A sub-sub-clustering.
#This is not benchmarked, but based on the same clustering logic
#that the more populated clusters make more sense,
#so use on your own risk :)

import os,sys,shutil
import glob
import math
from sklearn.cluster import DBSCAN
import numpy as np

def main(name, argv):
	if not len(argv) == 1:
		print_usage(name)
		return

	#Retrieving top resutls
	RMSD = float(argv[0])
	chain = 'B'
	top_folder = "Results/"
	files = [f for f in os.listdir() if '.pdb' in f]
	os.mkdir(top_folder)
	for f in files:
		shutil.copyfile(f, top_folder + f)
	#Clustering
	os.chdir(top_folder)
	big_clusters, num_labels = apply_DBSCAN(files, chain, RMSD)
	os.chdir('../')

def apply_DBSCAN(names, chain, threshold, cls_size_threshold=5):
	models = []
	for l in names:
		with open(l, 'r') as f:
			final = [[line[30:38], line[38:46], line[46:54]] for line in f if line[:4] == 'ATOM' and line[21] == chain and '  CA  ' in line]
		final = [[float(x) for x in line] for line in final]
		models.append(final)
	num = len(models)
	rmsd = []
	for i,m1 in enumerate(models):
		line = [0]*num
		for j,m2 in enumerate(models):
			for k in range(len(m1)):
				for l in [0,1,2]:
					line[j] += (m1[k][l]-m2[k][l])**2
			line[j] = math.sqrt(line[j]/len(m1)) 
		rmsd.append(line)
	rmsd = np.array(rmsd)
	clustering = DBSCAN(eps=threshold, metric='precomputed', min_samples=1).fit(rmsd)
	labels = clustering.labels_
	num_labels = max(labels) + 1
	cluster_size = [0]*num_labels
	cluster_sum = [0.0]*num_labels
	big_clusters = 0
	for i in range(num_labels):
		os.mkdir('tmp' + str(i + 1))
	for i,label in enumerate(labels):
		if not label == -1:
			cluster_size[label] += 1
			shutil.copyfile('../' + names[i], 'tmp' + str(label + 1) + '/' + names[i])
			os.remove(names[i])

	#Ranking the clusters, primarily by cluster size, then by average of final score
	clusters = zip(cluster_size, ['tmp' + str(i + 1) for i in range(num_labels)])
	sorted_clusters = sorted(clusters, reverse=True)
	for i,c in enumerate(sorted_clusters):
		os.rename(c[1], 'cluster' + str(i + 1))
		if c[0] >= cls_size_threshold:
			big_clusters += 1
			
	return big_clusters, num_labels

def print_usage(name):
	print("Usage : " + name + " <RMSD>")
if __name__ == "__main__":
	main(sys.argv[0], sys.argv[1:])
