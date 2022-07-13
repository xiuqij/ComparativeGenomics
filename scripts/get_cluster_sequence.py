import argparse

parser = argparse.ArgumentParser()
parser.add_argument('clusters', help='File with the identifiers')
parser.add_argument('proteomes', help='File with all protein sequences')

args = parser.parse_args()
def read_sequence(filename):
    '''Extract the header and sequence.
    Parameter filename: file with all protein sequences
    Yields: (header,sequence) for each entry
    '''
    header, seq = '', []

    # Locate the first sequence entry
    for line in filename:
        if line.startswith('>'):
            header, seq = line, []
            break
    # Read the rest of the sequences
    for line in filename:
        if line.startswith('>'):
            yield (header, ''.join(seq))
            header, seq = line, []
        else:
            seq.append(line)
    yield (header, ''.join(seq))

# Get sequence for each cluster
dict_identifiers = {}
with open(args.clusters,'r') as f:  # read the cluster identifiers into a dictionary
    clusterid = 1
    for line in f:
        identifiers = line.split('\t')
        dict_identifiers[clusterid] = identifiers
        clusterid += 1

    for clusterid, identifiers in dict_identifiers.items():
        outputname = "cluster_{id}.fasta".format(id = clusterid)    # create separate files for each cluster
        with open(outputname,'w') as f:
            for i in identifiers:
                proteomes = read_sequence(open(args.proteomes,'r'))
                for header,seq in proteomes:
                    if i.split('_')[1] in header and i.split('_')[2] in header:
                        f.write(header+seq)

                    



