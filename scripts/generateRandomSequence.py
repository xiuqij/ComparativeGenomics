import random
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('cluster', help = 'cluster file for determining how many sequences to generate')
parser.add_argument('output', help = 'path for output file')
args = parser.parse_args()

def random_aa_seq(length):
    '''
    create a random amino acid sequence at a given length.
    '''
    return ''.join(random.choice('ACDEFGHIKLMNPQRSTVWY') for i in range(length))

# Generate one 20-aa long sequence for each line in the cluster file
with open(args.output,'w') as output:
    with open(args.cluster,'r') as cluster:
        for line in cluster:    
            output.write(random_aa_seq(20)+'\n')