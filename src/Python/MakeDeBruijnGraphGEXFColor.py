#!/local/cluster/bin/python2.7
# !/usr/bin/env python
import sys, io, re, getopt, os, subprocess, math
#forget all this complicated stuff, we're just going to hard code it myself...
#sys.path.append('/raid2/labs/Koslicki_lab/koslickd/PythonScripts/pygraphml')
#from pygraphml.GraphMLParser import *
#from pygraphml.Graph import *
#from pygraphml.Node import *
#from pygraphml.Edge import *

def main(argv):
	input_file = ''
	kmer_size = 5
	initial_hash_size = '100M'
	num_threads = 5
	rev_comp_yes_no = 'no'
	lower_count = '0'
	out_file = ''
	r_val = '0'
	g_val = '0'
	b_val = '0'
	a_val = '.5'
	normalization = 'prob'
	try:
		opts, args = getopt.getopt(argv,"h:i:k:s:t:c:l:r:g:b:a:n:o:",["Help=","Inputfile=","KmerSize=","InitHash=", "NumThreads=","OutFile="])
	except getopt.GetoptError:
		print 'Unknown option, call using: ./MakeDeBruijnGraph.py -i <InputFile.fa> -k <KmerSize> -s <HashSize> -t <NumThreads> -c <ReverseComplementYesNo> -l <CountKmersAbove> -r -g -b -a <ColorSpecs0-255,0-1> -n <normalization, log, sqrt, linear, prob> -o <OutFile.gexf> '
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print './MakeDeBruijnGraph.py -i <InputFile.fa> -k <KmerSize> -s <HashSize> -t <NumThreads> -c <ReverseComplementYesNo> -l <CountKmersAbove> -r -g -b -a <ColorSpecs0-255,0-1> -n <normalization, log, sqrt, linear, prob> -o <OutFile.gexf>'
			sys.exit(2)
		elif opt in ("-i", "--InputFile"):
			input_file = arg
		elif opt in ("-k", "--KmerSize"):
			kmer_size = int(arg)
		elif opt in ("-s", "--HashSize"):
			initial_hash_size = arg
		elif opt in ("-t", "--NumThreads"):
			num_threads = int(arg)
		elif opt in ("-c", "--RevCompYesNo"):
			rev_comp_yes_no = arg
		elif opt in ("-l", "--LowerCount"):
			lower_count = arg
		elif opt in ("-o","--OutFile"):
			out_file = arg
		elif opt in ("-r","--Red"):
			r_val = arg
		elif opt in ("-g","--Green"):
			g_val = arg
		elif opt in ("-b","--Blue"):
			b_val = arg
		elif opt in ("-a","--Opacity"):
			a_val = arg
		elif opt in ("-n","--Normalization"):
			normalization = arg


	#Count the kmers in the file
	print('Counting ' + str(kmer_size) + 'mers using jellyfish to form vertices')
	if rev_comp_yes_no == 'no':
		cmd = './jellyfish count -m ' + str(kmer_size) + ' -L ' + lower_count + ' -s ' + initial_hash_size + ' -t '+ str(num_threads) +' ' + str(input_file) + ' -o ' + input_file + '-' + str(kmer_size) + 'mers.jf'
	else:
		cmd = './jellyfish count -m ' + str(kmer_size) + ' -L ' + lower_count + ' -C -s ' + initial_hash_size + ' -t '+ str(num_threads) +' ' + str(input_file) + ' -o ' + input_file + '-' + str(kmer_size) + 'mers.jf'
	result = subprocess.check_output(cmd, shell = True)

	#Dumpt the kmers to a file
	print('Dumping ' + str(kmer_size) + 'mers to file')
	cmd = './jellyfish dump ' + input_file + '-' + str(kmer_size) + 'mers.jf -o ' + input_file + '-' + str(kmer_size) + 'mers.fa'
	result = subprocess.check_output(cmd, shell = True)

	print('Counting ' + str(kmer_size+1) + 'mers using jellyfish to form edges')
	if rev_comp_yes_no == 'no':
		cmd = './jellyfish count -m ' + str(kmer_size+1) + ' -L ' + lower_count + ' -s ' + initial_hash_size + ' -t '+ str(num_threads) +' ' + str(input_file) + ' -o ' + input_file + '-' + str(kmer_size+1) + 'mers.jf'
	else:
		cmd = './jellyfish count -m ' + str(kmer_size+1) + ' -L ' + lower_count + ' -C -s ' + initial_hash_size + ' -t '+ str(num_threads) +' ' + str(input_file) + ' -o ' + input_file + '-' + str(kmer_size+1) + 'mers.jf'
	result = subprocess.check_output(cmd, shell = True)

	print('Dumping ' + str(kmer_size+1) + 'mers to file')
	cmd = './jellyfish dump ' + input_file + '-' + str(kmer_size+1) + 'mers.jf -o ' + input_file + '-' + str(kmer_size+1) + 'mers.fa'
	result = subprocess.check_output(cmd, shell = True)
	
	my_re = re.compile(r">") #Get rid of header file
	
	print('Reading in vertices')
	vertices = list()
	vertex_sizes = list()
	file_handle = open(input_file + '-' + str(kmer_size) + 'mers.fa','r')
	for line in file_handle:
		line_striped=line.strip()
		if not(my_re.search(line_striped)):
			vertices.append(line_striped)
		else:
			vertex_sizes.append(int(line_striped[1:]))
	file_handle.close()

	#Normalize the vertex weights
	total = sum(vertex_sizes)
	for index in range(0,len(vertex_sizes)):
		if normalization == 'prob':
			vertex_sizes[index] = float(vertex_sizes[index]) / float(total) 
		elif normalization == 'linear':
			vertex_sizes[index] = float(vertex_sizes[index])
		elif normalization == 'log':
			vertex_sizes[index] = math.log(float(vertex_sizes[index]))
		elif normalization == 'sqrt':
			vertex_sizes[index] = math.sqrt(float(vertex_sizes[index]))
		else:
			print('Unknown normalization choice, should be one of: prob, linear, log, sqrt')
	
	print('Reading in edges')
	edges = list()
	edge_sizes = list()
	file_handle = open(input_file + '-' + str(kmer_size + 1) + 'mers.fa','r')
	for line in file_handle:
		line_striped=line.strip()
		if not(my_re.search(line_striped)):
			edges.append(line_striped)
		else:
			edge_sizes.append(int(line_striped[1:]))
	file_handle.close()
	
	#Normalize edge weights
	total = sum(edge_sizes)
	for index in range(0,len(edge_sizes)):
		if normalization == 'prob':
			edge_sizes[index] = float(edge_sizes[index]) / float(total) 
		elif normalization == 'linear':
			edge_sizes[index] = float(edge_sizes[index])
		elif normalization == 'log':
			edge_sizes[index] = math.log(float(edge_sizes[index]))
		elif normalization == 'sqrt':
			edge_sizes[index] = math.sqrt(float(edge_sizes[index]))
		else:
			print('Unknown normalization choice, should be one of: prob, linear, log, sqrt')
	
	print('Forming graph')
	print('Writing vertices')
	file_handle = open(out_file,'w')
	file_handle.write('<?xml version="1.0" encoding="UTF-8"?>\n <gexf xmlns="http://www.gexf.net/1.2draft" xmlns:viz="http://www.gexf.net/1.1draft/viz" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.gexf.net/1.2draft http://www.gexf.net/1.2draft/gexf.xsd" version="1.2">\n \t <graph mode="static" defaultedgetype="directed">\n \t \t <nodes>\n')
	i=0
	for vertex in vertices:
		file_handle.write('\t \t \t <node id="%s" label="%s">\n' %(vertex,vertex))
		file_handle.write('\t \t \t \t <viz:size value="%s" />\n' % str(vertex_sizes[i]))
		file_handle.write('\t \t \t \t <viz:color r="%s" g="%s" b="%s" a="%s"/>\n' % (r_val,g_val,b_val,a_val))
		file_handle.write('\t \t \t </node>\n')
		i = i + 1
	file_handle.write('\t \t </nodes>\n')
	
	print('Writing edges')
	file_handle.write('\t \t <edges>\n')
	i=0
	for edge in edges:
		file_handle.write('\t \t \t <edge id="%s" source="%s" target="%s" >\n' %(edge,edge[0:len(edge)-1],edge[1:]))
		file_handle.write('\t \t \t \t <viz:size value="%s" />\n' % str(edge_sizes[i]))
		file_handle.write('\t \t \t \t <viz:color r="%s" g="%s" b="%s" a="%s"/>\n' % (r_val,g_val,b_val,a_val))
		file_handle.write('\t \t \t </edge>\n')
		i = i + 1
	file_handle.write('\t \t </edges>\n \t </graph>\n </gexf>')
	file_handle.close()
	
	print('Cleaning up')
	os.remove(input_file + '-' + str(kmer_size) + 'mers.jf')
	os.remove(input_file + '-' + str(kmer_size) + 'mers.fa')
	os.remove(input_file + '-' + str(kmer_size+1) + 'mers.jf')
	os.remove(input_file + '-' + str(kmer_size+1) + 'mers.fa')
        
    
if __name__ == "__main__":
	main(sys.argv[1:])