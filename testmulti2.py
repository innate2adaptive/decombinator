from __future__ import division
from Bio import SeqIO
from time import time, clock
from itertools import izip
from multiprocessing import Process, Queue
import sys
import argparse
import gzip
import os
import Levenshtein as lev
import collections as col
import multiprocessing as mp


# http://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file

fn = '2test'

def worker(arg, q):
    '''stupidly simulates long running process'''
    start = time.clock()
    s = 'this is a test'
    txt = s
    for i in xrange(200000):
        txt += s 
    done = time.clock() - start
    with open(fn, 'rb') as f:
        size = len(f.read())
    res = 'Process' + str(arg), str(size), done
    q.put(res)
    return res
  
def prnt(x1,x2,x3, q):
  #print str(len(x1.seq)), str(len(x2.seq)), str(len(x3.seq))
  outstr = str(x1.id) + "_" + str(x2.id) + "_" + str(x3.id)
  q.put(outstr)
  return outstr
  #count += 1
    

def listener(q):
    '''listens for messages on the q, writes to file. '''

    f = open(fn, 'wb') 
    while 1:
        m = q.get()
        if m == 'kill':
            f.write('killed')
            break
        f.write(str(m) + '\n')
        f.flush()
    f.close()

def main():
    #must use Manager queue here, or will not work
    manager = mp.Manager()
    q = manager.Queue()    
    pool = mp.Pool(4)

    #put listener to work first
    watcher = pool.apply_async(listener, (q,))

    #fire off workers
    jobs = []
    for i in range(80):
        #job = pool.apply_async(worker, (i, q))
        job = [pool.apply(prnt, args=(r1,r2,r3, q)) for r1,r2,r3 in izip(fq1, fq2, fq3)]
        jobs.append(job)
    return jobs
    # collect results from the workers through the pool result queue
    for job in jobs: 
        job.get()

    #now we are done, kill the listener
    q.put('kill')
    pool.close()


fq1 = SeqIO.parse(gzip.open("bigtestR1.fq.gz"), "fastq")
fq2 = SeqIO.parse(gzip.open("bigtestI1.fq.gz"), "fastq")
fq3 = SeqIO.parse(gzip.open("bigtestR2.fq.gz"), "fastq")



t0 = time() # Begin timer
  
main()


timed = time() - t0
took = round(timed,2)
#print count, 'reads processed from', rd1file, 'and', fq2file, 'and output into', outfq #FIX
if took < 60:
  print '\t\t\t\t\t\t\t\t\tTook', took, 'seconds to demultiplex samples'
else:
  print '\t\t\t\t\t\t\t\t\tTook', round((timed/60),2), 'minutes to jimmy indexes and hexamers around'









