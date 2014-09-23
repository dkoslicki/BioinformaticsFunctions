#!/local/cluster/bin/python2.7

import sys, io, re, getopt, os, subprocess, time

#This script will build a contig when given a seed k-mer

def main(argv):
	input_kmer=''
	counts_file=''
	lower_threshold = 0
	ave_crit = 'no'
	try:
		opts, args = getopt.getopt(argv,"h:k:c:L:A:m:",["Help=","Inputfile=","KmerSize=","InitHash=", "NumThreads=","OutFile="])
	except getopt.GetoptError:
		print 'Unknown option, call using: ./BuildContig.py -k <InputKmer> -c <CountsFile.jf> -L <LowerThreshold> -A <AverageCriterionYesNo> -m <MinLength>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print './BuildContig.py -k <InputKmer> -c <CountsFile.jf> -L <LowerThreshold> -A <AverageCriterionYesNo> -m <MinLength>'
			sys.exit(2)
		elif opt in ("-k", "--InputKmer"):
			input_kmer = arg
		elif opt in ("-c", "--CountsFile"):
			counts_file = arg
		elif opt in ("-L", "--LowerThreshold"):
			lower_threshold = int(arg)
		elif opt in ("-A","--AverageCriterionYesNo"):
			ave_crit = arg
		elif opt in ("-m", "--MinLength"):
			min_length = int(arg)

	#Always make sure that the input kmer is the same length as the counts file
	seq = input_kmer
	kmer_size = len(input_kmer)
	seen_kmers = dict()
	seen_kmers[input_kmer]=1

	#Get initial count
	cmd = './jellyfish query ' + counts_file +' '+ seq + ' | cut -d\' \' -f 2 '
	count_temp = int(subprocess.check_output(cmd, shell = True))
	counts = list()
	counts.append(count_temp)

	#Start the loop out with true conditions
	count = lower_threshold + 1
	i = 1
	#Extend to the right
	#print('Extending to the right')
	while count > lower_threshold:
		count = 0
		letter_to_add = ''
		for letter in ['A', 'C', 'T', 'G']:
			cmd = './jellyfish query ' + counts_file +' '+ seq[i:i + kmer_size - 1] + letter + ' | cut -d\' \' -f 2 '
			count_temp = int(subprocess.check_output(cmd, shell = True))
			#print(seq[i:i + kmer_size - 1] + letter)
			#See if this new count is the max or not
			if count_temp > count:
				count = count_temp
				letter_to_add = letter
		
		#If I've already seen the kmer, then just terminate. I.e. don't even try to deal with loops
		if seen_kmers.has_key(seq[i:i + kmer_size - 1] + letter):
			#print('Already seen kmer')
			break
		else:
			seen_kmers[seq[i:i + kmer_size - 1] + letter] = 1
		
		#Also, if the counts drop to half the average so far, then break out (too low of coverage)
		if ave_crit == 'yes' and count < (.5*sum(counts))/len(counts): #Make sure to do the float converseion
			break
		if count > 0:
			counts.append(count)
		seq += letter_to_add
		i = i + 1
	
	#Extend to the left
	#print('Extending to the left')
	count = lower_threshold + 1
	while count > lower_threshold:
		count = 0
		letter_to_add = ''
		for letter in ['A', 'C', 'G', 'T']:
			cmd = './jellyfish query ' + counts_file +' '+ letter + seq[0:kmer_size] + ' | cut -d\' \' -f 2 '
			count_temp = int(subprocess.check_output(cmd, shell = True))
			#print(letter + seq[0:kmer_size])
			#See if this new count is the max or not
			if count_temp > count:
				count = count_temp
				letter_to_add = letter

		#If I've already seen the kmer, then just terminate. I.e. don't even try to deal with loops
		if seen_kmers.has_key(letter + seq[0:kmer_size]):
			#print('Already seen kmer')
			break
		else:
			seen_kmers[letter + seq[0:kmer_size]] = 1

		#Also, if the counts drop to half the average so far, then break out (too low of coverage)
		if ave_crit == 'yes' and count < (.5*sum(counts))/len(counts): #Make sure to do the float converseion
			break
		if count > 0:
			counts[:0] = [count]
		seq = letter_to_add + seq

	if len(seq) >= min_length:
		print('>' + str(counts))
		print(seq)



if __name__ == "__main__":
	main(sys.argv[1:])