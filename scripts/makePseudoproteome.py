import argparse
import pandas as pd 

parser = argparse.ArgumentParser()

parser.add_argument('gene_order', help = 'Gene order file')
parser.add_argument('seq', help = 'random amino acid sequences')
parser.add_argument('output', help = 'path to output pseudoproteome file')
args = parser.parse_args()

# Read the random sequences and order them into a list
with open(args.seq,'r') as seq: 
    seq_list = []
    for line in seq:
        seq_list.append(line.rstrip())

# Read the gene order file and order them into a list
with open(args.gene_order,'r') as order:
    order_list = []
    for line in order:
        order_list.append(int(line.rstrip()))   # change the type for gene order to integer for sorting 

# Make a dataframe to assign sequences to gene order entries
data = {'Gene order':order_list, 'Sequence':seq_list}
df = pd.DataFrame(data)
# Sort the dataframe by gene order (numerical order)
df_sorted = df.sort_values(by = ['Gene order'], ignore_index = True)

# Write the sequences to the pseudoproteome file according to the sorted order
with open(args.output,'w') as output:
    output.write('>'+args.output+'\n')  # genome header
    for seq in df_sorted['Sequence']:
        output.write(seq)
    output.write('\n')  # Line break for combining multiple files later