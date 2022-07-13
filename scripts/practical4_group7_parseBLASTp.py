import argparse
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument('input_blast', help='The output file from BLASTp')
parser.add_argument('output', help='Output file with the best hits')
args = parser.parse_args()

# Read BLAST output 
with open(args.output,'w') as f:
    with open(args.input_blast, 'r') as blast:  # Read BLAST output line by line
        ref_orfs=[] # keep track of exsiting proteins
        for line in blast:
            record = line.split('\t')
            ref_orf = str(record[0])
            target_orf = str(record[1])
            
            # To extract only one hit for each orf 
            if ref_orf in ref_orfs:    #skip if the header already exists
                continue
            ref_orfs.append(ref_orf)
            
            f.write('Reference_'+ref_orf+'\t'+'Target_'+target_orf+'\n')