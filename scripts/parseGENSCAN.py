import argparse

parser = argparse.ArgumentParser()
parser.add_argument('genscan', help='The GENSCAN file to be processed')
parser.add_argument('output_cds', help='File to write the nucleotide sequences')
parser.add_argument('output_peptide', help='File to write the amino acid sequences')
args = parser.parse_args()

def read_sequence(filename):
    '''Extract the header and sequence of predicted peptide and coding sequences from the GENSCAN output file.
    Parameter filename: file with GENSCAN output
    Type filename: str

    Yields: (header,sequence) for each entry
    '''
    header, seq = None, []

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

# Reads input file and extract peptide and coding sequence
cds = open(args.output_cds,'w')
peptide = open(args.output_peptide,'w')

with open(args.genscan,'r') as f:
    for header, seq in read_sequence(f):
        if 'peptide' in header:
            peptide.write(header + seq)
        if 'CDS' in header:
            cds.write(header + seq) 

cds.close()
peptide.close()