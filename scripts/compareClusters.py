import argparse
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

parser = argparse.ArgumentParser()
parser.add_argument('cluster1', help='cluster file (BLASTp BBH)')
parser.add_argument('cluster2', help='cluster file (InParanoid)')
args = parser.parse_args()

# Define lists for storing cluster entries
cluster1 = []
cluster2 = []

# Read cluster file and store to corresponding list
with open(args.cluster1,'r') as f1:
    for line in f1:
        cluster1.append(line.rstrip())
    
with open(args.cluster2,'r') as f2:
    for line in f2:
        cluster2.append(line.rstrip())

# Compare the clusters and make Venn diagram
venn2([set(cluster1), set(cluster2)], ('Clusters by BLAST BBH', 'Clusters by InParanoid'))
plt.savefig('compare_cluster.png')
plt.show()