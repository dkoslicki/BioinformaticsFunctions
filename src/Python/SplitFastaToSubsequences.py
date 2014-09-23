#!/usr/bin/env python2.7
# !/local/cluster/bin/python2.7
#Untested
import sys, getopt
#sys.path.append("/local/cluster/biopython-1.58")
from Bio import SeqIO
import os
import re
def main(argv):
    input_file = ''
    out_seq_length = ''
    overlap_size = ''
    out_file = ''
    try:
        opts, args = getopt.getopt(argv,"hi:l:v:o:")
    except getopt.GetoptError:
        print 'Unknown option, call using: ./SplitFastaToSubsequences.py -i <InputFile.fa> -l <SequenceLength> -v <OverlapSize> -o <OutFile.fa>'
        sys.exit(2)
    if len(opts) <4:
        #print(opts)
        print('Too few options, call using: ./SplitFastaToSubsequences.py -i <InputFile.fa> -l <SequenceLength> -v <OverlapSize> -o <OutFile.fa>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'Call using: ./SplitFastaToSubsequences.py -i <InputFile.fa> -l <SequenceLength> -v <OverlapSize> -o <OutFile.fa>'
            sys.exit(2)
        elif opt in ("-i"):
            input_file = arg
        elif opt in ("-l"):
            out_seq_length = int(arg)
        elif opt in ("-v"):
            overlap_size = int(arg)
        elif opt in ("-o"):
            out_file = arg
    #In case the files are located outside of the directory where this script exits, let's get just the file name parts
    print('test')
    drive, input_file_path = os.path.splitdrive(input_file)
    input_file_path, input_file = os.path.split(input_file_path)
    drive, out_file_path = os.path.splitdrive(out_file)
    out_file_path, out_file = os.path.split(out_file_path)

    input_file_handle = open(os.path.join(input_file_path,input_file),'r')
    input_seqs = list(SeqIO.parse(input_file_handle,'fasta'))
    split_seqs = list()
    for seq in input_seqs:
        for subseq_index in range(0, len(seq) - out_seq_length,out_seq_length - overlap_size):
            #I can't seem to get the overlap right, so let's just make sure we don't go past the length of the sequence, perhaps it's right now...
            #if subseq_index + out_seq_length <= len(seq):
            split_seqs.append(seq[subseq_index:subseq_index + out_seq_length])

    out_file_handle = open(os.path.join(out_file_path,out_file),'w')
    SeqIO.write(split_seqs,out_file_handle,'fasta')

if __name__ == "__main__":
	main(sys.argv[1:])