import argparse

parser = argparse.ArgumentParser()
parser.add_argument('blastp_1', help='BLASTp output ordered by reference ORFs')
parser.add_argument('blastp_2', help='BLASTp output ordered by target ORFs')
parser.add_argument('output', help='Output file with the best hits')
args = parser.parse_args()

# Read the blast outputs and sort the BHs as dictionaries
with open(args.blastp_1, 'r') as blastp1:  # Read BLAST output line by line
    ref_orfs=[] # keep track of exsiting proteins
    ref_tar_BH = {} # store the BHs with genes in reference proteome as key
    for line in blastp1:
        record = line.split('\t')
        ref_orf = str(record[0])
        target_orf = str(record[1])           
        # To extract only one hit for each orf 
        if ref_orf in ref_orfs:    #skip if the header already exists
            continue
        ref_orfs.append(ref_orf)    # add to exsiting proteins
        ref_tar_BH[ref_orf] = target_orf    # add as BH pair

with open(args.blastp_2, 'r') as blastp2:
    tar_orfs=[] # keep track of exsiting proteins
    tar_ref_BH = {} # store the BHs with genes in target proteome as key
    for line in blastp2:
        record = line.split('\t')
        target_orf = str(record[0])
        ref_orf = str(record[1])

        if target_orf in tar_orfs:  #skip if the header already exists
            continue
        tar_orfs.append(target_orf) # add to exsiting proteins
        tar_ref_BH[target_orf] = ref_orf    # add as BH pair

# Find BBH and write to output
with open(args.output, 'w') as output:
    for ref, tar in ref_tar_BH.items():
        if tar_ref_BH[tar] == ref:  # checkes if this BH is also the BH in the target proteome
            output.write('Reference_' + ref + '\t' + 'Target_' + tar + '\n')