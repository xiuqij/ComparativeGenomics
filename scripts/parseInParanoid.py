import argparse
import itertools

parser = argparse.ArgumentParser()
parser.add_argument('InParanoid', help='InParanoid output table')
parser.add_argument('output', help='File to write clustered output')
args = parser.parse_args()

# Group proteins with the same group-id
with open(args.InParanoid,'r') as ip:
    proteins = {}
    for line in ip:
        record = line.split('\t')
        groupid = record[0]
        protein = record[4]
        if groupid not in proteins:
            proteins[groupid] = []  # create a new item if the groupid was not in the current dictionary
        proteins[groupid].append(protein.rstrip())
# Get all ortholog pairs and write to output
with open(args.output,'w') as output:
    for groupid, protein_list in proteins.items():
        combinations = itertools.combinations(protein_list,2)
        for i in combinations:
            if i[0].split('_')[0] != i[1].split('_')[0]:    # If the two proteins belong to different organism
                output.write('Reference_' + i[0] +'\t'+ 'Target_' +i[1]+'\n')