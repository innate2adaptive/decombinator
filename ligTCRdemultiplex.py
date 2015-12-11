# ligTCRdemultiplex.py v1.0
# James M. Heather, December 2015, UCL

##################
### BACKGROUND ###
##################

# A derivate of DualIndexDemultiplexing.py/FuzzyDID, fuzzily demultiplexes ligation TCR sequencing protocol FASTQ data 
  # i.e. allows a specified number of mismatches in the index sequence
# Takes all three reads and simultaneously demultiplexes and formats read 1 for vDCR.py analysis
  # i.e. adds the first 30 nucleotides of R2 to R1
  # This will contain the random barcode followed by any V(D)J information downstream, allowing collapsing
# Demultiplexes using a combination of indexes from third read and another index nested in read 1

##################
###### INPUT #####
##################

# Requires the command line input of at least 3file names, giving the three Illumina reads
  # Files may be uncompressed, or gzipped (and be named accordingly, e.g. File.fastq.gz)
# A fourth optional comma delimited file detailing sample index specifics is strongly recommended, allowing production of correctly named files
  # File must give the following details, one sample (or index combination) per line, with no empty lines:
    # Sample name, SP1/R1 index (I), SP2/R2 index (L):
      # e.g.: P005v1,1,11
# e.g. run: python ligTCRdemultiplex.py -r1 read1.fastq -r2 read2.fastq -i1 indexread1.fastq -ix indexes.ndx

# Other optional flags:
  
  # -s/--summary: Output a summary file containing details of the run into a 'Logs' directory. Default = True
  
  # -a/--outputall: Output the results of all possible index combinations currently used in protocol
    # e.g. Useful in finding potential cross-contaminating or incorrectly indexed samples
    # NB - This option can be run even if an index list is provided (although only those provided by the index list will be named)
  
  # -t/--threshold: Specifies the threshold by which indexes can be clustered by fuzzy string matching, allowing for sequencing errors
    # Default = 2. Setting to zero turns off fuzzy matching, i.e. only allowing exact string matching
  
  # -z/--gzip: Specify whether or not to gzip compress the output, demultiplexed FASTQs. Default = True
  
  # -c/--count: Specify whether or not to show the running line count, every 100,000 reads. 
    # Helps in monitoring progress of large batches. Default = True.

# To see all options, run: python ligTCRdemultiplex.py -h

##################
##### OUTPUT #####  
##################
    
# A fastq file will be produced for each sample listed in the index file, in the modified format, containing all reads that matched
# So we go from:        R1 - [6s|X1|----VDJ-----]
#                       R2 - [X2]
#                       R3 - [8s|N1|8s|N2|2s|-----5'UTR-----]
# To: ========>         out- [8s|N1|8s|N2|2s|X1|X2|----VDJ-----]         
# Where X = hexamer index sequence, N = random barcode sequence, and Ys = spacer sequences of length Y
  # The 8s sequences can be used downstream to identify the presence and location of random barcode N sequences
  # 2s is also kept to allow for the possibility finding N sequences produced from slipped reads

##################
#### PACKAGES ####  
##################

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

##########################################################
############# READ IN COMMAND LINE ARGUMENTS #############
##########################################################

def args():
  """args(): Obtains command line arguments which dictate the script's behaviour"""

  # Help flag
  parser = argparse.ArgumentParser(
      description='Script to demultiplex FASTQ data produced using the ligation TCRseq protocol')
  # Add arguments
  parser.add_argument(
      '-r1', '--read1', type=str, help='Read 1 FASTQ file', required=True)
  parser.add_argument(
      '-r2', '--read2', type=str, help='Read 1 FASTQ file', required=True)
  parser.add_argument(
      '-i1', '--index1', type=str, help='Index read FASTQ file', required=True)
  parser.add_argument(
      '-ix', '--indexlist', type=str, help='File containing sample/index table', required=False)
  parser.add_argument(
      '-s', '--summary', type=bool, help='Output summary data (True/False)', required=False, default=True)
  parser.add_argument(
      '-a', '--outputall', type=bool, help='Output all possible index combinations (True/False)', required=False, default=False)
  parser.add_argument(
      '-t', '--threshold', type=int, help='Edit distance allowed for fuzzy-string index matching (default=2)', required=False, default=2)
  parser.add_argument(
      '-z', '--gzip', type=bool, help='Gzip compress output FASTQ files (True/False)', required=False, default=True)
  parser.add_argument(
      '-c', '--count', type=bool, help='Show the count (True/False)', required=False, default=True)
  

  return parser.parse_args()

inputargs = vars(args())

############# FIX - remember to add a check so that if output all = true and indexlist given, it will produce named output for the given ones

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

outputfiles = ["Undetermined.fq"]

# If given an indexlist, use that to generate named output files
if inputargs['indexlist']:
  indexes = list(open(inputargs['indexlist'], "rU"))
else:
  print "No index list file provided."
  sys.exit()

XXdict = {}

for x in indexes:
  
  if x == "\n":
    print "Empty line detected in index file, presumed end of file."
    break

  elements = x.strip("\n").split(",")

  sample = elements[0]
  
  open(sample + ".fq", "w").close()
  
  compound_index = X1dict[elements[1]] + X2dict[elements[2]] 
  
  XXdict[compound_index] = open(sample + ".fq", "a")
  
  outputfiles.append(sample + ".fq")
  

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

print "Reading input files..."

# Open read files
if inputargs['read1'].endswith('.gz'):
  fq1 = SeqIO.parse(gzip.open(inputargs['read1']), "fastq")
else:
  fq1 = SeqIO.parse(open(inputargs['read1']), "fastq")
  
if inputargs['index1'].endswith('.gz'):
  fq2 = SeqIO.parse(gzip.open(inputargs['index1']), "fastq")
else:
  fq2 = SeqIO.parse(open(inputargs['index1']), "fastq")


if inputargs['read2'].endswith('.gz'):
  fq3 = SeqIO.parse(gzip.open(inputargs['read2']), "fastq")
else:
  fq3 = SeqIO.parse(open(inputargs['read2']), "fastq")

print "Demultiplexing data..."

for record1, record2, record3 in izip(fq1, fq2, fq3):
  
  count += 1  

  if count % 100000 == 0 and inputargs['count'] == True:
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
      
      if lev.distance(ndx, seqX) <= inputargs['threshold']:
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

# Gzip compress output

print "Compressing demultiplexed files..."

if inputargs['gzip'] == True:
  for f in outputfiles:  
    
    with open(f) as infile, gzip.open(f + '.gz', 'wb') as outfile:
        outfile.writelines(infile)
    os.unlink(f)

##########################################################
#################### STDOUT STATISTICS ###################
##########################################################

timed = time() - t0
took = round(timed,2)
#print count, 'reads processed from', rd1file, 'and', fq2file, 'and output into', outfq #FIX
if took < 60:
  print '\t\t\t\t\t\t\t\t\tTook', took, 'seconds to demultiplex samples'
else:
  print '\t\t\t\t\t\t\t\t\tTook', round((timed/60),2), 'minutes to jimmy indexes and hexamers around'


print count, "reads processed"
print dmpd_count, "reads demultiplexed"
print fuzzy_count, "reads demultiplexed using fuzzy index matching"

if clash_count > 0:
  print clash_count, "reads had fuzzy index clashes (i.e. could have assigned to >1 index) and were discarded"


# Write out list of fuzzy matched sequences, so can fish out later if needed

if inputargs['threshold'] > 0:
  
  print "\nOutputting list of reads demultiplexed using fuzzy index matching"

  fuzzout = open("FuzzyMatchedIDs.txt", "w")
  for f in fuzzies:
    print >> fuzzout, f

  fuzzout.close()
else:
  
  print "Fuzzy index matching not used (threshold = 0)"

