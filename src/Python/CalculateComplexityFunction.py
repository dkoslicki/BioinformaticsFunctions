#!/local/cluster/bin/python2.7
import sys, getopt, subprocess, os, re, resource, numpy
sys.path.append("/local/cluster/biopython-1.58")
from Bio import SeqIO

#This script will calculate the complexity function (number of unique kmers) for a given fasta file

def main(argv):
	input_file = ''
	split_yes_no = 'no'
	kmer_size = ''
	prefix_file_name = ''
	out_file = ''
	initial_hash_size = '100M'
	max_kmer = 101
	min_kmer = 10
	num_jellyfish_threads = 5
	already_run = 'no'
	try:
		opts, args = getopt.getopt(argv,"h:i:s:j:m:M:o:",["InputFile=","SplitYesNo="])
	except getopt.GetoptError:
		print 'Unknown option, call using: ./CalculateComplexityFunction.py -i <InputFile.fa> -s <HashSize> -j <JellyfishThreads> -m <min_kmer> -M <max_kmer> -o <Outfile.csv>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print './CalculateComplexityFunction.py -i <InputFile.fa> -s <HashSize> -m <min_kmer> -j <JellyfishThreads> -m <MaxKmer> -o <Outfile.csv>'
			sys.exit(2)
		elif opt in ("-i", "--InputFile"):
			input_file = arg
		elif opt in("-k", "--KmerSize"):
			kmer_size = int(arg)
		elif opt in ("-s", "--InitHash"):
			initial_hash_size = arg
		elif opt in ("-j", "--JellyfishThreads"):
			num_jellyfish_threads = int(arg)
		elif opt in ("-M", "--MaxKmer"):
			max_kmer = int(arg)
		elif opt in ("-o","--OutFile"):
			out_file = arg
		elif opt in ("-m","--MinKmer"):
			min_kmer = int(arg)
			
	
	out_vect = numpy.empty([1, max_kmer],dtype=numpy.int)
	for kmer_size in range(min_kmer, max_kmer + 1):
		print('Counting ' + str(kmer_size) + 'mers')
		#Calvin's code is  much faster for smaller kmer sizes
		if kmer_size >= 16:
			cmd = './jellyfish count -m ' + str(kmer_size) + ' -s ' + initial_hash_size + ' -t '+ str(num_jellyfish_threads) +' ' + input_file + ' --out-counter 1 -o ' + input_file +'-' + str(kmer_size) + 'mers.jf'
			result = subprocess.check_output(cmd, shell = True)
		
			print('Counting number of unique ' + str(kmer_size) + 'mers')
			#It's way faster to have jellyfish output the histogram and then sum up the counts, as this will give the number of unique kmers
			cmd = './jellyfish histo ' + input_file +'-' + str(kmer_size) + "mers.jf | cut -d' ' -f2 | awk '{s+=$1}END{print s}'"
			num_kmers = int(subprocess.check_output(cmd, shell = True))
			os.remove(input_file +'-' + str(kmer_size) + 'mers.jf')
		else:
			cmd = 'kmer_total_count -i ' + input_file + ' -k ' + str(kmer_size) + ' -n -l | wc -l'
			num_kmers = int(subprocess.check_output(cmd, shell = True))

		out_vect[0, kmer_size - 1] = num_kmers
		
	numpy.savetxt(out_file, out_vect, delimiter=',',fmt='%d')

if __name__ == "__main__":
	main(sys.argv[1:])