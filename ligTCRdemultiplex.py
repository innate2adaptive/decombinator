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
  
  # -s/--supresssummary: Supress the production of a summary file containing details of the run into a 'Logs' directory. Default = False
  
  # -a/--outputall: Output the results of all possible index combinations currently used in protocol
    # e.g. Useful in finding potential cross-contaminating or incorrectly indexed samples
    # NB - This option can be run even if an index list is provided (although only those provided by the index list will be named)
  
  # -t/--threshold: Specifies the threshold by which indexes can be clustered by fuzzy string matching, allowing for sequencing errors
    # Default = 2. Setting to zero turns off fuzzy matching, i.e. only allowing exact string matching
  
  # -z/--dontgzip: Suppress the automatic compression of output demultiplexed FASTQ files with gzip. Default = False
  
  # -c/--dontcount: Suppress the whether or not to show the running line count, every 100,000 reads. 
    # Helps in monitoring progress of large batches. Default = False.

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
from itertools import izip
import time
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
      '-t', '--threshold', type=int, help='Edit distance allowed for fuzzy-string index matching (default=2)', required=False, default=2)
  parser.add_argument(
      '-s', '--suppresssummary', type=bool, help='Output summary data (True/False)', required=False, default=False)
  parser.add_argument(
      '-a', '--outputall', type=bool, help='Output all possible index combinations (True/False)', required=False, default=False)
  parser.add_argument(
      '-z', '--dontgzip', type=bool, help='Stop the output FASTQ files automatically being compressed with gzip (True/False)', required=False, default=False)
  parser.add_argument(
      '-c', '--dontcount', type=bool, help='Show the count (True/False)', required=False, default=False)
  parser.add_argument(
      '-fl', '--fuzzylist', type=bool, help='Output a list of thosereads that demultiplexed using fuzzy index matching (True/False)', required=False, default=False)  
  parser.add_argument(
      '-ex', '--extension', type=str, help='Specify the file extension of the output FASTQ files. Default = \"fq\"', required=False, default="fq")
  return parser.parse_args()

inputargs = vars(args())

if inputargs['outputall'] == False and not inputargs['indexlist']:
  print "No indexing file provided, and output all option not enabled; one (or both) is required."
  sys.exit()

##########################################################
############ CREATE DICTIONARIES FOR INDEXES #############
##########################################################

# SP1 index = R1 (our own, RC1 proximal index)
X1dict = {"1":"ATCACG", "2":"CGATGT", "3":"TTAGGC", "4":"TGACCA", "5":"ACAGTG", "6":"GCCAAT", "11":"GGCTAC", "12":"CTTGTA", "13":"TAGACT"}
# NB: One index removed due to similarity to others, but exists in earlier datasets: "14":"ACACGG" 


# 'SP2' index = R2 (index read, comes first in rearranged sequence)
X2dict = {"1":"CGTGAT", "2":"ACATCG", "3":"GCCTAA", "4":"TGGTCA", "5":"CACTGT", "6":"ATTGGC", "7":"GATCTG", "8":"TCAAGT", "9":"CTGATC", "10":"AAGCTA", "11":"GTAGCC", "12":"TACAAG", "13":"TTGACT", "14":"GGAACT"}

##########################################################
########### GENERATE SAMPLE-NAMED OUTPUT FILES ###########
##########################################################


suffix = "." + inputargs['extension']

failed = open("Undetermined" + suffix, "w")

outputreads = col.Counter()
outputreads["Undetermined"] = 0

usedindexes = col.defaultdict(list)       # This keeps a track of all files that have been generated to house demultiplexed reads

XXdict = {}


# If given an indexlist, use that to generate named output files
if inputargs['indexlist']:
  indexes = list(open(inputargs['indexlist'], "rU"))

  for x in indexes:
    
    if x == "\n":
      print "Empty line detected in index file, presumed end of file."
      break

    elements = x.strip("\n").split(",")
    sample = elements[0]
    
    open(sample + suffix, "w").close()
    
    compound_index = X1dict[elements[1]] + X2dict[elements[2]] 
    XXdict[compound_index] = open(sample + suffix, "a")
    
    outputreads[sample] = 0
    usedindexes[sample] = compound_index


# If the outputall option is chosen, output all possible index combinations that exist in the data
  # Note that if an indexlist is provided, those names are still used in the appropriate output files
  # Also note that while all combinations are looked for, those which remain unused at the end will be deleted
  
if inputargs['outputall'] == True:
  # generate all possible index combinations, and then check if they have been generated yet (via an index file)
    # only make those that haven't
  allXcombs = [x + "-" + y for x in X1dict.keys() for y in X2dict.keys()]
  
  for x in allXcombs:
    compound_index = X1dict[x.split("-")[0]] + X2dict[x.split("-")[1]] 
    
    if compound_index not in usedindexes.values():
      XXdict[compound_index] = open("Indexes_" + x + suffix, "a")
      outputreads["Indexes_" + x] = 0
      usedindexes["Indexes_" + x] = compound_index
      

count = 0
dmpd_count = 0          # number successfully demultiplexed 
fuzzy_count = 0		# number of sequences that were demultiplexed using non-exact index matches
clash_count = 0		# number of fuzzy ID clashes

fuzzies = []		# list to record IDs matched using fuzzy indexes

t0 = time.time() # Begin timer
  
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

  if count % 100000 == 0 and inputargs['dontcount'] == False:
    print '\t read', count
  
### NB For non-standard Illumina encoded fastqs, might need to change which fields are carried into fq_* vars
  
  fq_id = record1.id

  # N relates to barcode random nucleotides, X denotes index bases
  
  ### FORMATTING OUTPUT READ ###
  
  Nseq = record3.format("fastq").split('\n')[1][0:30]
  X1seq = record1.format("fastq").split('\n')[1][6:12]
  X2seq = str(record2.seq)
  readseq = str(record1.seq)[12:]

  Nqual = record3.format("fastq").split('\n')[3][0:30]
  X1qual = record1.format("fastq").split('\n')[3][6:12]
  X2qual = record2.format("fastq").split('\n')[3]
  readqual = record1.format("fastq").split('\n')[3][12:]
    
  fq_seq = Nseq + X1seq + X2seq + readseq
  fq_qual = Nqual + X1qual + X2qual + readqual
  
  new_record = str("@" + fq_id + "\n" + fq_seq + "\n+\n" + fq_qual + "\n")  
  
  ### DEMULTIPLEXING ###
  
  seqX = X1seq + X2seq 
  
  if seqX in XXdict:
    # Exact index matches
          
    XXdict[seqX].write(new_record)
    
    dmpd_count += 1
    outputreads[str(XXdict[seqX]).split("\'")[1][:len(suffix)]] += 1
    
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
      outputreads[str(XXdict[matches[0]]).split("\'")[1][:len(suffix)]] += 1
      
    else:
      
      if len(matches) > 1:
	clash_count += 1
	
      failed.write(new_record)
      outputreads['Undetermined'] += 1
  
for x in XXdict.values():
  x.close()

failed.close()

# If output all is allowed, delete all unused index combinations
if inputargs['outputall'] == True:
  for f in outputreads.keys():
    if outputreads[f] == 0:
      os.remove(f + suffix)
      del outputreads[f]
      del usedindexes[f]

# Gzip compress output

if inputargs['dontgzip'] == False:
  print "Compressing demultiplexed files..."
  
  for f in outputreads.keys():  
    
    with open(f + suffix) as infile, gzip.open(f + suffix + '.gz', 'wb') as outfile:
        outfile.writelines(infile)
    os.unlink(f + suffix)

#################################################
################## STATISTICS ###################
#################################################

timed = time.time() - t0
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

# Write data to summary file
if inputargs['suppresssummary'] == False:
  
  # Check for directory and make summary file
  if not os.path.exists('Logs'):
    os.makedirs('Logs')
  date = time.strftime("%Y_%m_%d")
  summaryfile = open("Logs/" + date + "_Demultiplexing_Summary.csv", "w")
  
  # Generate string to write to summary file
  
  summstr = "Property,Value\nDirectory," + os.getcwd() + "\nDateFinished," + date + "\nTimeFinished," + time.strftime("%H:%M:%S") + "\nTimeTaken(Seconds)," + str(round(timed,2)) + "\n"
  
  for s in ['read1', 'read2', 'index1', 'indexlist', 'extension', 'threshold', 'outputall', 'dontgzip', 'fuzzylist']:
    summstr = summstr + s + "," + str(inputargs[s]) + "\n"
  
  summstr = summstr + "NumberReadsInput," + str(count) + "\nNumberReadsDemultiplexed," + str(dmpd_count) + "\nNumberFuzzyDemultiplexed," + str(fuzzy_count) + "\nNumberIndexClash," + str(clash_count) + "\n\nOutputFile,IndexUsed\n" 
  
  
  # Write out number of reads in and details of each individual output file
  for x in sorted(usedindexes.keys()):
    summstr = summstr + x + "," + usedindexes[x] + "\n"
  
  summstr = summstr + "\nOutputFile,NumberReads\n"
  
  for x in sorted(outputreads.keys()):
    summstr = summstr + x + "," + str(outputreads[x]) + "\n"
    
  print >> summaryfile, summstr 
  
  summaryfile.close()
  
# Write out list of fuzzy matched sequences, so can fish out later if needed
if inputargs['threshold'] > 0 and inputargs['fuzzylist'] == True:
  
  print "\nOutputting list of reads demultiplexed using fuzzy index matching"
  
  # Check for directory and make summary file
  if not os.path.exists('Logs'):
    os.makedirs('Logs')
  date = time.strftime("%Y_%m_%d")
  
  fuzzout = open("Logs/" + date + "_FuzzyMatchedIDs.txt", "w")
  for f in fuzzies:
    print >> fuzzout, f

  fuzzout.close()

