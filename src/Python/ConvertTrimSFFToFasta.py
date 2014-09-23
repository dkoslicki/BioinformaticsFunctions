#!/local/cluster/bin/python2.7
import sys
sys.path.append("/local/cluster/biopython-1.58")
from Bio import SeqIO
import os
import re
all_file_names = os.listdir('/raid2/labs/Koslicki_lab/shirlemo/SRR05')
for file_name in all_file_names:
	file_name_parts = re.split('\.',file_name)
	base_name = file_name_parts[0]
	extension = file_name_parts[1]
	if extension == "sff":
		file_handle = open(file_name,'r')
		seqs = list(SeqIO.parse(file_handle,"sff-trim"))
		file_handle.close()
		out_handle = open(base_name + ".fasta","w")
		SeqIO.write(seqs,out_handle,"fasta")
		out_handle.close()

#Next, from the command line: ls *.fasta | head -n 1 > AllFiles.fa
#Then: ls *.fasta | xargs -I{} cat {} AllFiles.fa >> AllFiles.fa