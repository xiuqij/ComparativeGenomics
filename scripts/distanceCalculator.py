import argparse
import math
import itertools

parser = argparse.ArgumentParser()
parser.add_argument('input', help='The file with GC content of each genome')
#parser.add_argument('output_matrix', help='Output file with the distance matrix')
args = parser.parse_args()

def calc_dist(gc1,gc2):
    '''Calculate the distance between the two given GC values.'''
    dist = math.sqrt((gc1-gc2) ** 2)
    return dist

# Read the input file with GC content of each genome and extract the values into a dictionary
with open(args.input,'r') as input:
    gc_dict = {}    # with genome as key and GC content as the value
    for line in input:
        genome = line.split(' ')[0]
        gc_content = float(line.split(' ')[1].rstrip())
        gc_dict[genome] = gc_content

# Calculate the distance for each pair
distance_dict = {}  # with the genome pair as key, and their distance as value
combs = itertools.combinations(list(gc_dict.keys()),2)  # generate combinations of different pairs
for i in combs:
    distance = calc_dist(gc_dict[i[0]],gc_dict[i[1]])
    distance_dict[i] = distance

# print result
for pair, dist in distance_dict.items():
    print("{}: {}\n".format(pair,dist))