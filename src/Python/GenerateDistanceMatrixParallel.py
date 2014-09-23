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


kmer_size = 10
num_threads = 40

Alph = ['A', 'C', 'G', 'T']
E = list(product(Alph, repeat=kmer_size))
E = [list(elem) for elem in E]
E = ["".join(elem) for elem in E]
kmers_list = E
print(len(kmers_list))
D = numpy.zeros([len(kmers_list),len(kmers_list)],dtype = numpy.int8)

#Compute it in parallel
pool = Pool(processes = num_threads)
rows = pool.map(compute_row_star,izip(kmers_list, repeat(kmers_list), repeat(kmer_size)))
out_file_h5 = h5py.File('D' + str(kmer_size) + '.h5', 'w')
out_file_h5.create_dataset("data", data = D, compression="gzip", compression_opts=9)
