#!/local/cluster/bin/python2.7

import sys, io, re, getopt, os, subprocess, h5py, numpy
numpy.set_printoptions(threshold=numpy.nan)
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
	input_fileA = 'uninitialized'
	input_fileB = 'uninitialized'
	output_file = 'uninitialized'
	kmer_size = 5
	num_threads = 5
	rev_comp_yes_no = 'no'
	lower_cutoff = 1
	jellyfish_threads = 1
	hash_size = '100M'
	try:
		opts, args = getopt.getopt(argv,"h:A:B:t:L:o:",["Help=","Inputfile=","KmerSize=","InitHash=", "NumThreads=","OutFile="])
	except getopt.GetoptError:
		print 'Unknown option, call using: ./GenerateCountsAndReducedDistanceMatrix.py -A <InputCounts1.fa> -B <InputCounts2.fa> -t <NumThreads> -L <LowerCountCutoff> -o <OutputFile.h5>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print './GenerateCountsAndReducedDistanceMatrix.py -A <InputCounts1.fa> -B <InputCounts2.fa> -t <NumThreads> -L <LowerCountCutoff> -o <OutputFile.h5>'
			sys.exit(2)
		elif opt in ("-A", "--InputFile1"):
			input_fileA = arg
		elif opt in ("-B", "--InputFile2"):
			input_fileB = arg
		elif opt in ("-o", "--OutputFile"):
			output_file = arg
		elif opt in ("-t", "--NumThreads"):
			num_threads = int(arg)
		elif opt in ("-k", "--KmerSize"):
			kmer_size = int(arg)
		elif opt in ("-c", "--RevCompYesNo"):
			rev_comp_yes_no = arg
		elif opt in ("-L", "--LowerCutoff"):
			lower_cutoff = int(arg)
		elif opt in ("-j", "--JellyfishThreads"):
			jellyfish_threads = int(arg)
		elif opt in("-s", "HashSize"):
			hash_size = arg


	#This assumes we've already done the kmer counts, now we just read them in. Don't use a lower cutoff when making the counts, I can handle that right here in this code

	#First, read in the kmers (should be produced via jellyfish, and then jellyfish dump), it MUST have the format '>count\n kmer\n'
	second_flag = 0
	for input_file in [input_fileA, input_fileB]:
		print('Reading in kmers and counts')
		vertices = list()
		vertex_sizes = list()
		file_handle = open(input_file,'r')
		i = 1
		for line in file_handle:
			if i % 2 == 1:
				line_stripped_split = line.strip()
				#print(line_stripped_split[1:])
				temp_count = int(line_stripped_split[1:])
				if temp_count >= lower_cutoff: #only add it if the associated count is above the threshold
					vertex_sizes.append(temp_count)
					count_flag = 1
				else:
					count_flag = 0
			elif count_flag:
				line_stripped_split = line.strip()
				#print(line_stripped_split)
				vertices.append(line_stripped_split)
			i = i + 1

		file_handle.close()
		if second_flag:
			B_kmers = vertices
			B_sizes = vertex_sizes
		else:
			A_kmers = vertices
			A_sizes = vertex_sizes
			second_flag = 1
	
	kmer_size = len(A_kmers[0])

	#Compute the union
	print('Forming the union')
	my_union = dict()
	for kmer in A_kmers:
		my_union[kmer] = 1
	for kmer in B_kmers:
		my_union[kmer] = 1

	kmers_list = my_union.keys()

	print('Forming the count vectors')
	A_dict = dict(zip(A_kmers,A_sizes))
	B_dict = dict(zip(B_kmers,B_sizes))
	
	A_counts = numpy.zeros([1,len(kmers_list)],dtype = numpy.int64)
	B_counts = numpy.zeros([1,len(kmers_list)],dtype = numpy.int64)

	print('Forming first vector')
	for i in range(len(kmers_list)):
		if A_dict.has_key(kmers_list[i]):
			A_counts[0,i] = A_dict[kmers_list[i]]

	print('Forming second vector')
	for i in range(len(kmers_list)):
		if B_dict.has_key(kmers_list[i]):
			B_counts[0,i] = B_dict[kmers_list[i]]
	
	print('Saving results')
	numpy.savetxt(input_fileA + '-' + input_fileB + '-' + str(kmer_size) + 'mer-counts.csv', A_counts, fmt = '%u', delimiter = ',')
	numpy.savetxt(input_fileB + '-' + input_fileA + '-' + str(kmer_size) + 'mer-counts.csv', B_counts, fmt = '%u', delimiter = ',')

	print('Saving Kmers')
	file_handle = open(input_fileA + '-' + input_fileB + '-' + str(kmer_size) + 'mers.txt','w')
	for i in range(len(kmers_list)):
		file_handle.write('%s\n' % kmers_list[i])
	file_handle.close()
		

	print('Number of Kmers: ' + str(len(kmers_list)))

	print('Forming distance matrix')

	#Initialize the alphabet, kmer size, and distance matrix
	Alph = ['A', 'C', 'G', 'T']
	D = numpy.zeros([len(kmers_list),len(kmers_list)],dtype = numpy.int8)

	#Compute it in parallel
	pool = Pool(processes = num_threads)
	rows = pool.map(compute_row_star,izip(kmers_list, repeat(kmers_list), repeat(kmer_size)))

	print('Combining results')
	for index in range(len(kmers_list)):
		D[index,:] = rows[index]

	print('Saving results')
	out_file_h5 = h5py.File(output_file, 'w')
	out_file_h5.create_dataset("data", data = D, compression="gzip", compression_opts=9)
#	out_file_h5.create_dataset("data", data = D)
	numpy.savetxt(output_file + '.txt', D, fmt='%d', delimiter=',')




	

if __name__ == "__main__":
	main(sys.argv[1:])
