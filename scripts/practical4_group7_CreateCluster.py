import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input1', help='File with best hit alignments')
parser.add_argument('input2', help='File with best hit alignments')
parser.add_argument('input3', help='File with best hit alignments')
parser.add_argument('output', help='File to write clustered output')
args = parser.parse_args()

clusters = {}   # with reference protein as key, and target list as the value
with open(args.input1,'r') as f1:
    for line in f1:
        ref = line.split('\t')[0]
        target = line.split('\t')[1]
        if ref not in clusters:
            clusters[ref] = []
        clusters[ref].append(target.rstrip())

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

with open(args.output,'w') as cluster_file:
    for ref, target in clusters.items():
        if len(target) == 3:    # Extract only the items with matches from all three target proteomes 
            target_list = '\t'.join(target)
            cluster_file.write(ref+'\t'+target_list+'\n')