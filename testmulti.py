from __future__ import division
from Bio import SeqIO
from time import time, clock
from itertools import izip
from multiprocessing import Pool
import sys
import argparse
import gzip
import os
import Levenshtein as lev
import collections as col
import multiprocessing

def worker():
    """worker function"""
    print 'Worker'
    return


def prnt(x1,x2,x3):
  #print str(len(x1.seq)), str(len(x2.seq)), str(len(x3.seq))
  outstr = str(x1.id) + "_" + str(x2.id) + "_" + str(x3.id)
  print >>outfile, outstr
  #count += 1
  
  
  
  #print str(len(record1.seq)), str(len(record2.seq)), str(len(record3.seq))

# Open read files
fq1 = SeqIO.parse(gzip.open("testR1.fq.gz"), "fastq")
fq2 = SeqIO.parse(gzip.open("testI1.fq.gz"), "fastq")
fq3 = SeqIO.parse(gzip.open("testR2.fq.gz"), "fastq")

outfile = open("testout.fq", "w")

p = Pool(4)
#y = p.map(prnt, (for record1, record2, record3 in izip(fq1, fq2, fq3)))

#for record1, record2, record3 in izip(fq1, fq2, fq3):
  #print str(len(record1.seq)), str(len(record2.seq)), str(len(record3.seq))

##jobs = []
##for i in range(5):
    ##p = multiprocessing.Process(target=worker)
    ##jobs.append(p)
    ##p.start()
    
    
#p.map(prnt, args=(r1,r2,r3)) for x in range(1,7)
global count
count = 0

t = [p.apply(prnt, args=(r1,r2,r3,)) for r1,r2,r3 in izip(fq1, fq2, fq3)]
p.close()
outfile.close()

print str(len(t))
print str(count)


