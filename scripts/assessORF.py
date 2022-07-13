import argparse

parser = argparse.ArgumentParser()
parser.add_argument('orf_ref', help='BLASTp output of predicted ORFs against reference proteome')
parser.add_argument('ref_orf', help='BLASTp output of reference proteome against predicted ORFs')
args = parser.parse_args()

TP_count = 0
FP_count = 0
FN_count = 0
significant_cutoff = 0.001

# Read the BLASTp output ordered by predicted ORFs
with open(args.orf_ref, 'r') as blast1:
    orfs=[] # keep track of exsiting entries
    for line in blast1:
        record = line.split('\t')

        orf = str(record[0])
        if orf in orfs:    # only check the best hit for each ORF
            continue
        orfs.append(orf)

        evalue = float(record[10])  # extract e-value
        if evalue <= significant_cutoff:  # if the hit is a significant match
            TP_count += 1
        else:   # if the hit is not a significant match
            FP_count += 1

# Read the BLASTp output ordered by the entries in the reference proteome
with open(args.ref_orf, 'r') as blast2:
    refs=[] # keep track of exsiting proteins
    for line in blast2:
        record = line.split('\t')

        ref = str(record[0])
        if ref in refs:    # only check the best hit for each protein
            continue
        refs.append(ref)

        evalue = float(record[10])
        if evalue > significant_cutoff:  # if this protein does not have a significant match in the predicted ORFs
            FN_count += 1  

# Calculate precision and recall
precision = TP_count/(TP_count+FP_count)
recall = TP_count/(TP_count+FN_count)

# print result
print("Precision: {:.2f}; Recall: {:.2f}.".format(precision,recall))