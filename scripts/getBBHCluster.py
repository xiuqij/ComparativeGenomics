import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input1', help='File with best hit alignments')
parser.add_argument('input2', help='File with best hit alignments')
parser.add_argument('input3', help='File with best hit alignments')
parser.add_argument('output', help='File to write clustered output')
args = parser.parse_args()

clusters = {}   # key: reference protein; value: a list of corresponding target proteins
with open(args.input1,'r') as f1:
    for line in f1:
        ref = line.split('\t')[0]   # read the reference protein
        target = line.split('\t')[1]    # read the target protein
        if ref not in clusters:
            clusters[ref] = []  # create a new item if the reference protein was not in the current dictionary
        clusters[ref].append(target.rstrip())   # append the target protein to the list

# do the same for the rest of the file
with open(args.input2,'r') as f2:
    for line in f2:
        ref = line.split('\t')[0]
        target = line.split('\t')[1]
        if ref not in clusters:
            clusters[ref] = []
        clusters[ref].append(target.rstrip())

with open(args.input3,'r') as f3:
    for line in f3:
        ref = line.split('\t')[0]
        target = line.split('\t')[1]
        if ref not in clusters:
            clusters[ref] = []
        clusters[ref].append(target.rstrip())
# write to the cluster file
with open(args.output,'w') as cluster_file:
    for ref, target in clusters.items(): 
        target_list = '\t'.join(target)
        cluster_file.write(ref+'\t'+target_list+'\n')