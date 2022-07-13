import argparse
parser = argparse.ArgumentParser()
parser.add_argument('cluster', help = 'input cluster file')
parser.add_argument('pfa', help = 'original pfa file')
parser.add_argument('output', help = 'path for output file')
args = parser.parse_args()

# Give index to the proteins according to their order in pfa file
with open(args.pfa,'r') as pfa: # Give an index to each entry in the pfa file
    orf_order = {}
    index = 0   # starts with 0
    for line in pfa:
        if line.startswith('>'):    # extract header
            orf = line.split('_')[1]    # extract orf id
            orf_order[index] = orf.rstrip()
            index += 1
    key_list = list(orf_order.keys())   # store keys and values to list 
    value_list = list(orf_order.values())

# Read cluster file and extract the ORF IDs in their order of appearance 
with open(args.cluster,'r') as clusters:   
    orf_list = []
    # Read the cluster file and extract records according to the genome number
    if '17' in args.pfa:
        for line in clusters:
            record = line.split('\t')
            orf_list.append(record[0].split('_')[2].rstrip())   # extract ORF name
    if '20' in args.pfa:        
        for line in clusters:
            record = line.split('\t')
            orf_list.append(record[1].split('_')[2].rstrip())
    if '23' in args.pfa:        
        for line in clusters:
            record = line.split('\t')
            orf_list.append(record[2].split('_')[2].rstrip())
    if '26' in args.pfa:        
        for line in clusters:
            record = line.split('\t')
            orf_list.append(record[3].split('_')[2].rstrip())

# Find the corresponding index for the ORFs and write to output
with open(args.output,'w') as output:   
    for header in orf_list:
        index = value_list.index(header)
        order_number = key_list[index]  # Retrive the order in the pfa file
        output.write(str(order_number) + '\n')