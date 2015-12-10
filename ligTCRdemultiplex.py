# FuzzyDid.py v1.0
# James M. Heather, May 2015, UCL

##################
### BACKGROUND ###
##################

# A derivate of DualIndexDemultiplexing.py, that allows fuzzy index matching
  # i.e. allows a specified number of mismatches in the index sequence
# Takes all three reads and simultaneously demultiplexes and formats read 1 for vDCR.py analysis
  # i.e. the first 12 nucleotides of the read make up the random barcode, with any V(D)J information downstream
# Demultiplexes using a combination of indexes from third read and another index nested in read 1

##################
###### INPUT #####
##################

# Takes command line input of 4 file names:
  # First 3 are (Illumina encoded) fastq files, relating to reads 1, 2 (index read) and 3 respectively
  # Fourth is a comma delimited file detailing sample index specifics, one per line:
    # Sample name, SP1/R1 index (I), SP2/R2 index (L):
    # eg: P005v1,1,11
# Run: python DualIndexDemultiplexing.py read1.fastq read2.fastq read3.fastq indexes.ndx

##################
##### OUTPUT #####  
##################
    
# A fastq file will be produced for each sample listed in the index file, in the modified format, containing all reads that matched
# So we go from:        R1 - [N1|X1|----VDJ-----]
#                       R2 - [X2]
#                       R3 - [N2|-----5'UTR-----]
# To: ========>         out- [N2|N1|X1|X2|----VDJ-----]         

##################
#### PACKAGES ####  
##################

from __future__ import division
from Bio import SeqIO
from time import time, clock
from itertools import izip
import sys
import os
import Levenshtein as lev

####################
#### THRESHOLD! ####  
####################

# This sets the edit distance threshold allowed for fuzzy indexing

threshold = 2

# Users can also choose whether or not to output a list of those reads that were demultiplexed using fuzzy matching

fuzz_ids = True

##########################################################
############# READ IN COMMAND LINE ARGUMENTS #############
##########################################################

filename = ""

if (len(sys.argv) <> 5):
  print "Please supply the three read fastq files and sample/index file (e.g. python DualIndexDemultiplexing.py read1.fastq read2.fastq read3.fastq indexes.ndx)"
  sys.exit()
else:
  rd1file = str(sys.argv[1])
  rd2file = str(sys.argv[2])
  rd3file = str(sys.argv[3])
  indexfile = str(sys.argv[4])

fq1 = SeqIO.parse(open(rd1file), "fastq")
fq2 = SeqIO.parse(open(rd2file), "fastq")
fq3 = SeqIO.parse(open(rd3file), "fastq")


### Create dictionaries of the two kinds of indices

##########################################################
############ CREATE DICTIONARIES FOR INDEXES #############
##########################################################

# SP1 index = R1 (our own, RC1 proximal index)
X1dict = {"1":"ATCACG", "2":"CGATGT", "3":"TTAGGC", "4":"TGACCA", "5":"ACAGTG", "6":"GCCAAT", "11":"GGCTAC", "12":"CTTGTA", "13":"TAGACT", "14":"ACACGG"}

# 'SP2' index = R2 (index read, comes first in rearranged sequence)
X2dict = {"1":"CGTGAT", "2":"ACATCG", "3":"GCCTAA", "4":"TGGTCA", "5":"CACTGT", "6":"ATTGGC", "7":"GATCTG", "8":"TCAAGT", "9":"CTGATC", "10":"AAGCTA", "11":"GTAGCC", "12":"TACAAG", "13":"TTGACT", "14":"GGAACT"}

##########################################################
########### GENERATE SAMPLE-NAMED OUTPUT FILES ###########
##########################################################

failed = open("Undetermined.fq", "w")

indexes = list(open(indexfile, "rU"))

XXdict = {}

for x in indexes:

  elements = x.strip("\n").split(",")

  sample = elements[0]
  
  open(sample + ".fq", "w").close()
  
  compound_index = X1dict[elements[1]] + X2dict[elements[2]] 
  
  XXdict[compound_index] = open(sample + ".fq", "a")
  
  ## ADD IF INDEX ERROR STOP -> ENSURE NO EMPTY LINE AT END OF INDEX FILE

count = 0
dmpd_count = 0          # number successfully demultiplexed 
fuzzy_count = 0		# number of sequences that were demultiplexed using non-exact index matches
clash_count = 0		# number of fuzzy ID clashes

fuzzies = []		# list to record IDs matched using fuzzy indexes

t0 = time() # Begin timer
  

##########################################################
########### LOOP THROUGH ALL READ FILES IN SYNC ##########
######## PROCESS INTO CORRECT FORMAT & DEMULTIPLEX #######
##########################################################

for record1, record2, record3 in izip(fq1, fq2, fq3):
  
  count += 1  

  if count % 100000 == 0:
    print '\t read', count
  
### NB For non-standard Illumina encoded fastqs, might need to change which fields are carried into fq_* vars
  
  fq_id = record1.id

  # N relates to barcode random nucleotides, X denotes index bases
  
  ### FORMATTING OUTPUT READ ###
  
  N1seq = record1.format("fastq").split('\n')[1][0:6]
  N2seq = record3.format("fastq").split('\n')[1][0:6]
  X1seq = record1.format("fastq").split('\n')[1][6:12]
  X2seq = str(record2.seq)
  readseq = str(record1.seq)[12:]

  N1qual = record1.format("fastq").split('\n')[3][0:6]
  N2qual = record3.format("fastq").split('\n')[3][0:6]
  N2qual = record3.format("fastq").split('\n')[3][0:6]
  X1qual = record1.format("fastq").split('\n')[3][6:12]
  X2qual = record2.format("fastq").split('\n')[3]
  readqual = record1.format("fastq").split('\n')[3][12:]
    
  fq_seq = N2seq + N1seq + X1seq + X2seq + readseq
  fq_qual = N2qual + N1qual + X1qual + X2qual + readqual
  
  new_record = str("@" + fq_id + "\n" + fq_seq + "\n+\n" + fq_qual + "\n")  
  
  ### DEMULTIPLEXING ###
  
  seqX = X1seq + X2seq 
  
  if seqX in XXdict:
    # Exact index matches
          
    XXdict[seqX].write(new_record)
    
    dmpd_count += 1
    
  else:
    # Otherwise allow fuzzy matching
    
    matches = []
    
    for ndx in XXdict.keys():
      
      if lev.distance(ndx, seqX) <= threshold:
        matches.append(ndx)
    
    if len(matches) == 1:
      # Only allow fuzzy match if there is one candidate match within threshold
      XXdict[matches[0]].write(new_record)
      dmpd_count += 1
      fuzzy_count += 1
      fuzzies.append(fq_id)
      
    else:
      
      if len(matches) > 1:
	clash_count += 1
	
      failed.write(new_record)
  
for x in XXdict.values():
  x.close()

failed.close()

##########################################################
#################### STDOUT STATISTICS ###################
##########################################################

timed = time() - t0

#print count, 'reads processed from', rd1file, 'and', fq2file, 'and output into', outfq #FIX
print '\t\t\t\t\t\t\t\t\tTook', round(timed,2), 'seconds to jimmy indexes and hexamers around'
print count, "reads processed"
print dmpd_count, "reads demultiplexed"
print fuzzy_count, "reads demultiplexed using fuzzy index matching"

if clash_count > 0:
  print clash_count, "reads had fuzzy index clashes (i.e. could have assigned to >1 index) and were discarded"


# Write out list of fuzzy matched sequences, so can fish out later if needed

if fuzz_ids == True:
  
  print "\nOutputting list of reads demultiplexed using fuzzy index matching"

  fuzzout = open("FuzzyMatchedIDs.txt", "w")
  for f in fuzzies:
    print >> fuzzout, f

  fuzzout.close()

