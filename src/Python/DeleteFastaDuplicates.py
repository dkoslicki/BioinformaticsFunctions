#!/local/cluster/bin/python2.7

import sys, io, re, getopt, os, subprocess, time
sys.path.append("/local/cluster/biopython-1.58")
from Bio import SeqIO

#This script will uniquify a fasta file

def main(argv):
	input_file=''
	output_file=''
	try:
		opts, args = getopt.getopt(argv,"h:i:o:",["Help=","Inputfile=","KmerSize=","InitHash=", "NumThreads=","OutFile="])
	except getopt.GetoptError:
		print 'Unknown option, call using: ./DeleteFastaDuplicates.py -i <InputFile.fa> -o <Outfile.fa>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print './DeleteFastaDuplicates.py -i <InputFile.fa> -o <Outfile.fa>'
			sys.exit(2)
		elif opt in ("-i", "--InputFile"):
			input_file = arg
		elif opt in ("-o", "--OutputFile"):
			output_file = arg

	
	

	my_re = re.compile(r">") #Get rid of header file
	sequences = list()
	headers = list()
	file_handle = open(input_file,'r')	
	for line in file_handle:
		line_striped=line.strip()
		if not(my_re.search(line_striped)):
			sequences.append(line_striped)
		else:
			headers.append(line_striped)
	file_handle.close()

	sequence_to_header = dict()
	for i in range(len(sequences)):
		seq = sequences[i]
		header = headers[i]
		if not(sequence_to_header.has_key(seq)):
			sequence_to_header[seq]=header

	out_seqs = sequence_to_header.keys()
	out_headers = sequence_to_header.values()
	file_handle = open(output_file,'w')
	for i in range(len(out_seqs)):
		file_handle.write('%s\n' % out_headers[i]) 
		file_handle.write('%s\n' % out_seqs[i])
	file_handle.close()

if __name__ == "__main__":
	main(sys.argv[1:])