import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('input_genome', help='The Genome file to be processed')
parser.add_argument('input_prot', help='The Proteome file to be processed')
parser.add_argument('output_statistics', help='File to write the GC content, di-nucleotide frequencies and di-aminoacid frequencies')
args = parser.parse_args()

def check_seq(sequence):
	'''Checkes if the sequence contains undefined nucleotides'''
	undefined_nt_count = 0
	for i in sequence:
		if i not in 'ATGC':
			undefined_nt_count += 1
	if undefined_nt_count != 0:
		print("The sequence contains {} undefined nucleotides.".format(undefined_nt_count))
	else:
		print("The sequence does not contain undefined nucleotides.")

def calc_GC(sequence):
	'''Calculate the GC content of a given genome sequence'''
	A_count = sequence.count('A')
	T_count = sequence.count('T')
	G_count = sequence.count('G')
	C_count = sequence.count("C")
	gc_content = float(G_count+C_count)/float(A_count+T_count+G_count+C_count)
	return gc_content

def calc_dinucleotide(sequence):
	'''Calculate the dinucleotide counts of a given genome sequence'''
	dinucleotide_count = {}
	for i in range(len(sequence)-1):
		dinucleotide = sequence[i:i+2]
		if dinucleotide in dinucleotide_count:
			dinucleotide_count[dinucleotide] += 1
		else:
			dinucleotide_count[dinucleotide] = 1
	return dinucleotide_count

diamino_count = {}	# defined outside the function for combining the result from multiple sequence
def calc_diamino(sequence):
	'''Calculate the di- amino acid counts of a given proteome'''
	for i in range(len(sequence)-1):
		diamino = sequence[i:i+2]
		if diamino in diamino_count:
			diamino_count[diamino] += 1
		else:
			diamino_count[diamino] = 1

## Read the genome file and extract the nucleotide sequence
with open(args.input_genome, "r") as genome:
    genome_seq = ''
    for line in genome:
        if line.startswith('>'):
            continue
        else:
            genome_seq += line.rstrip()

# Calculate GC content and the dinucleotide frequencies
check_seq(genome_seq)
GC = calc_GC(genome_seq)
dint_count = calc_dinucleotide(genome_seq)

## Read the proteome file and extract protein sequences
with open(args.input_prot,'r') as proteome:
	prot_seqs = []
	seq = ''
	# Locate the first entry
	for line in proteome:
		if line.startswith('>'):
			break
    # Read the rest of the sequences
	for line in proteome:
		if line.startswith('>'):
			prot_seqs.append(seq)
			seq = '' # reset the current sequence
			continue
		else:
			seq += line.rstrip()
	prot_seqs.append(seq)	# for the last sequence

# Calculate diamino acid frequencies
for prot_seq in prot_seqs:
	calc_diamino(prot_seq)

# Total dinucleotide count and total diamino acid count
sum_dint = 0
sum_diamino = 0
for dint, count in dint_count.items():
	sum_dint += count
for diamino, count in diamino_count.items():
	sum_diamino += count

# Write output
with open(args.output_statistics,'w') as output:
	output.write("GC content	{:.2f}".format(GC) + '\n')
	output.write("Di-nucleotide frequency" + '\n')
	for dint, count in dint_count.items():
		output.write("{} = {}/{}".format(str(dint),count,sum_dint) + '\n')
	output.write("Di-amino acid frequency" + '\n')
	for diamino, count in diamino_count.items():
		output.write("{} = {}/{}".format(str(diamino),count,sum_diamino) + '\n')

# Plot the dinucleotide frequency
dints = list(dint_count.keys())
dint_count = list(dint_count.values())
plt.bar(dints,dint_count)
figname1 = str(args.input_genome) + '_dint.png'
plt.savefig(figname1)

# Plot the diamino acid frequency
diaminos = list(diamino_count.keys())
diamino_count = list(diamino_count.values())
plt.bar(diaminos,diamino_count)
figname2 = str(args.input_genome) + '_diamino.png'
plt.savefig(figname2)