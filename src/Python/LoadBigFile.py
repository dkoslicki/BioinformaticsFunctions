#!/local/cluster/bin/python2.7
import sys
import resource
sys.path.append("/local/cluster/biopython-1.58")
from Bio import SeqIO
import os
import re
file_handle = open('/data/temp/WGS_Database_NCBI_nt.fa','r')
print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000)
seqs = list(SeqIO.parse(file_handle,"fasta"))
print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000)
