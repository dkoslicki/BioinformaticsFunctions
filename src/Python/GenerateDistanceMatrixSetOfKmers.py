#!/local/cluster/bin/python2.7

import sys, io, re, getopt, os, subprocess, h5py, numpy
from itertools import *
from multiprocessing import Pool, freeze_support

def compute_row(kmer_i,kmers_list,kmer_size):
	Drow = numpy.zeros([1,len(kmers_list)],dtype = numpy.int8)
	i=0
	for kmer_j in kmers_list:
		for t in range(kmer_size + 1):
			if kmer_i[t:kmer_size] == kmer_j[0:kmer_size - t]:
				break
		Drow[0,i] = t
		i = i + 1
	return Drow

def compute_row_star(arg):
	return compute_row(*arg)


def main(argv):
	input_file = 'uninitialized'
	output_file = 'uninitialized'
	try:
		opts, args = getopt.getopt(argv,"h:i:t:o:",["Help=","Inputfile=","KmerSize=","InitHash=", "NumThreads=","OutFile="])
	except getopt.GetoptError:
		print 'Unknown option, call using: ./GenerateDistanceMatrixSetOfKmers.py -i <InputKmerFile.fa> -t <NumThreads> -o <OutputFile.h5>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print './GenerateDistanceMatrixSetOfKmers.py -i <InputKmerFile.fa> -t <NumThreads> -o <OutputFile.h5>'
			sys.exit(2)
		elif opt in ("-i", "--InputFile"):
			input_file = arg
		elif opt in ("-o", "--OutputFile"):
			output_file = arg
		elif opt in ("-t", "--NumThreads"):
			num_threads = int(arg)


	#First, read in the kmers (should be produced via jellyfish, and then jellyfish dump)
	my_re = re.compile(r">") #Get rid of header file


	print('Reading in Kmers')
	kmers_list = list()
	file_handle = open(input_file,'r')
	for line in file_handle:
		line_striped=line.strip()
		if not(my_re.search(line_striped)):
			kmers_list.append(line_striped)
	file_handle.close()

	print('Number of Kmers: ' + str(len(kmers_list)))

	print('Forming distance matrix')

	#Initialize the alphabet, kmer size, and distance matrix
	Alph = ['A', 'C', 'G', 'T']
	kmer_size = len(kmers_list[0])
	D = numpy.zeros([len(kmers_list),len(kmers_list)],dtype = numpy.int8)

	#Compute it in parallel
	pool = Pool(processes = num_threads)
	for kmer_i in kmers_list:
		rows = pool.map(compute_row_star,izip(kmers_list, repeat(kmers_list), repeat(kmer_size)))

#Serial
#	rows = list()
#	for index in range(len(kmers_list)):
#		rows.append(compute_row(kmers_list[index],kmers_list,kmer_size))

	
	for index in range(len(kmers_list)):
		D[index,:] = rows[index]


	out_file_h5 = h5py.File(output_file, 'w')
	out_file_h5.create_dataset("data", data = D, compression="gzip", compression_opts=9)
#	out_file_h5.create_dataset("data", data = D)



	

if __name__ == "__main__":
	main(sys.argv[1:])
