import argparse
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument('input_blast', help='The output file from BLAST')
parser.add_argument('input_db', help='The file containing all genome sequences')
parser.add_argument('output', help='Output file with the top scoring hits from each genome')
args = parser.parse_args()

# Create list for output sequences
output_sequences = []

# Read BLAST output and extract sequences
with open(args.input_blast, 'r') as blast:  # Read BLAST output line by line
    for line in blast:
        record = line.split('\t')
        genome_header = str(record[1])
        sbegin = int(record[8])
        send = int(record[9])
        
        # To extract only one hit for each genome 
        headers = []
        for i in output_sequences:
            headers.append(i.id)
        if genome_header in headers:    #skip if the header already exists
            continue

        # Extract the sequence according to the coordinates
        genome_sequences = SeqIO.parse(open(args.input_db), "fasta")
        for genome in genome_sequences:
            if genome_header == genome.id:  # find the matching entry
                if sbegin > send: # if the coordinates are in reverse order
                    sbegin, send = send, sbegin
                    sbegin -= 1
                    extracted_seq = genome.seq[sbegin:send]
                    extracted_seq = extracted_seq.reverse_complement()
                    output_sequences.append(SeqRecord(extracted_seq, id=genome.id, description=''))
                else: # if the coordinates are not in reverse order
                    sbegin -= 1
                    extracted_seq = genome.seq[sbegin:send]
                    output_sequences.append(SeqRecord(extracted_seq, id=genome.id, description=''))
            
# Write the sequences to a fasta file
SeqIO.write(output_sequences, open(args.output,'w'), "fasta")