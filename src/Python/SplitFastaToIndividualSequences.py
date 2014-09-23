#!/local/cluster/bin/python2.7
import sys, getopt
import resource
sys.path.append("/local/cluster/biopython-1.58")
from Bio import SeqIO
import os
import re

def main(argv):
	input_file = ''
	my_prefix = ''
	out_file = ''
	try:
		opts, args = getopt.getopt(argv,"h:i:p:",["InputFile=","Prefix="])
	except getopt.GetoptError:
		print 'Unknown option, call using: ./SplitFastaToIndividualSequences.py -i <InputFile.fa> -p <Prefix>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print './SplitFastaToIndividualSequences.py -i <InputFile.fa> -p <Prefix>'
			sys.exit(2)
		elif opt in ("-i", "--InputFile"):
			input_file = arg
		elif opt in("-p", "--Prefix"):
			my_prefix = str(arg)
	print(my_prefix)

	drive, input_path = os.path.splitdrive(input_file)
	input_path, garbage = os.path.split(input_path)
	file_handle = open(input_file,'r')
	print('Reading in file')
	seqs = list(SeqIO.parse(file_handle,'fasta'))
	file_handle.close()
	print('Writing files')
	i = 1
	for seq in seqs:
		file_handle = open(os.path.join(input_path,my_prefix + '_' + str(i) + '.fa'),'w')
		SeqIO.write(seq,file_handle,'fasta')
		file_handle.close()
		i = i + 1
            
            
if __name__ == "__main__":
	main(sys.argv[1:])