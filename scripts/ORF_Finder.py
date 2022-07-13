import argparse
import itertools
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument('input', help = 'input genome file (.fasta)')
parser.add_argument('minimum', type = int, default = 100, help = 'Minimum length cutoff of the ORFs')
parser.add_argument('output', help = 'output file with the ORF sequence')
args = parser.parse_args()

class ORF:
    '''An ORF object contains the sequence and other information on the identified ORFs.'''
    def __init__(self, genome, orf_id, coord, sequence, reverse):
        self.genome = genome    # input genome
        self.orf_id = orf_id    # an assigned identifier for the ORF
        self.coord = coord  # starting coordinate of the ORF
        self.sequence = sequence    # ORF sequence
        self.length = len(self.sequence)    # the length of the ORF
        self.reverse = reverse  # is the ORF on the reverse strand (True/False)

    def translate_sequence(self):
        '''Translates the ORF into protein sequence'''
        return self.sequence.translate()

    def create_header(self):
        '''Generates the header for the ORF when writing to output file'''
        if self.reverse:
            return ">{}_ORF{}_rev".format(self.genome,self.orf_id)
        else:
            return ">{}_ORF{}".format(self.genome,self.orf_id)

raw_orf_list = []   # used to store ORFs found on both strands
def search_orf(seq,reverse = False):
    '''
    Scan the input sequence in different reading frames and record all possible ORFs in a list.
    '''
    stop_codons = ['TAA','TAG','TGA']
    start_codon = 'ATG'
    frames = [0,1,2]    # read the sequence from these three frames
    orf_found = False  # tracks if an ORF is currently being read
    
    position = 0    # start position of the ORF
    orf_id = 1  # identifier for the ORF
    ORFsequence = ''    # stores the ORF sequence

    for i in frames:
        while i < len(seq)-2: 
            codon = seq[i:(i+3)]
            if codon == start_codon:    # if start codon is found, store the ORF that follows
                orf_found = True
                position = i+1  # record start position of the ORF
                ORFsequence += codon
                i += 3
            while orf_found == True:
                codon = seq[i:(i+3)]
                if codon not in stop_codons:
                    ORFsequence += codon
                    i += 3
                    if i > len(seq):    # if no stop codon is encountered until the end of the sequence, discard the current ORF
                        ORFsequence = ''    # reset ORF sequence
                        orf_found = False   # reset status
                        break
                else:   # if a stop codon is encountered, store the current ORF to a new ORF object
                    raw_orf_list.append(ORF(args.input,orf_id,position,ORFsequence,reverse))
                    ORFsequence = ''        # reset ORF sequence
                    orf_found = False   # reset status
                    orf_id += 1 # create new identifier
            i += 3

def find_overlap(orf_list):
    '''Finds overlapped genes in a given list of ORF objects, and returns the list containing the shorter ORFs with overlap.'''
    overlapped_genes = []
    comb = itertools.combinations(orf_list,2)   # create combinations of any two ORFs and check if they overlap
    for i in comb:
        orf1 = i[0]
        orf2 = i[1]
        if orf1.reverse == orf2.reverse:    # in the two ORFs are on the same strand
            if orf1.coord < orf2.coord: # compare starting position
                if orf1.coord + orf1.length > orf2.coord:   # if the two ORFs overlap
                    short_orf = orf1 if orf1.length < orf2.length else orf2 # take the shorter ORF (to be removed)
                    overlapped_genes.append(short_orf)
            else:
                if orf2.coord + orf2.length > orf1.coord:
                    short_orf = orf1 if orf1.length < orf2.length else orf2
                    overlapped_genes.append(short_orf)
    return overlapped_genes

def filter_orf(raw_orfs):
    '''Filters out: 1) the shorter ones in overlapped ORFs 2) ORFs that are shorter than the stated cutoff'''
    filtered_orf_list = []
    overlap_orfS = find_overlap(raw_orfs)
    for orf in raw_orfs:
        if (orf.length >= args.minimum) and (orf not in overlap_orfS):
            filtered_orf_list.append(orf)
    return filtered_orf_list

# Read genome file and get sequence
with open(args.input,'r') as genome:
    seq = ''
    for line in genome:
        if line.startswith('>'):
            continue
        else:
            seq += line.rstrip()
# Get input sequence on both strand
genome_seq = Seq(seq)   # forward strand
genome_seq_rc = genome_seq.reverse_complement() # reverse strand
# Get raw ORFs
search_orf(genome_seq)
search_orf(genome_seq_rc,reverse=True)
# Get filtered ORFs
filtered_orfs = filter_orf(raw_orf_list)
# Write the translated ORFs to output file
with open(args.output,'w') as output:
    for orf in filtered_orfs:
        header = orf.create_header()
        peptide = orf.translate_sequence()
        output.write(header+'\n'+str(peptide)+'\n')