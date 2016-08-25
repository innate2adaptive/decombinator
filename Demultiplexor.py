# Demultiplexor
# James M. Heather, August 2016, UCL
# https://innate2adaptive.github.io/Decombinator/

##################
### BACKGROUND ###
##################

# Takes all three reads and simultaneously demultiplexes and formats read 1 for vDCR.py analysis
  # i.e. adds the first 30 nucleotides of R2 to R1
  # This will contain the random barcode followed by any V(D)J information downstream, allowing collapsing
# Demultiplexes using a combination of indexes from third read and another index nested in read 1

# A derivate of DualIndexDemultiplexing.py/FuzzyDID, fuzzily demultiplexes ligation TCR sequencing protocol FASTQ data 
  # i.e. allows a specified number of mismatches in the index sequence
  
##################
###### INPUT #####
##################

# Requires the command line input of at least 3 file names, giving the three Illumina reads
  # Files may be uncompressed, or gzipped (and be named accordingly, e.g. File.fastq.gz)
# A fourth optional comma delimited file detailing sample index specifics is strongly recommended, allowing production of correctly named files
  # File must give the following details, one sample (or index combination) per line, with no empty lines:
    # Sample name, SP1/R1 index (I), SP2/R2 index (L):
      # e.g.: AlphaSample,1,11
    # If you include one and only one chain description (i.e. alpha, beta, gamma or delta) into your sample name, you need not set chain in Decombinator
    
# e.g. run: python Demultiplexor.py -r1 read1.fastq -r2 read2.fastq -i1 indexread1.fastq -ix indexes.ndx

# Other optional flags:
  
  # -s/--supresssummary: Supress the production of a summary file containing details of the run into a 'Logs' directory. 
  
  # -a/--outputall: Output the results of all possible index combinations currently used in protocol
    # e.g. Useful in finding potential cross-contaminating or incorrectly indexed samples
    # NB - This option can be run even if an index list is provided (although only those provided by the index list will be named)
  
  # -t/--threshold: Specifies the threshold by which indexes can be clustered by fuzzy string matching, allowing for sequencing errors
    # Default = 2. Setting to zero turns off fuzzy matching, i.e. only allowing exact string matching
  
  # -dz/--dontgzip: Suppress the automatic compression of output demultiplexed FASTQ files with gzip. 
    # Using this flag makes the script execute faster, but data will require more storage space. 
    
  # -dc/--dontcount: Suppress the whether or not to show the running line count, every 100,000 reads. 
    # Helps in monitoring progress of large batches. 
    
  # -fz/--fuzzylist: Output a list of FASTQ IDs of reads which are demultiplexed using fuzzy (i.e. non-exact) index matching, within the specified threshold.
    # Default = False, but can be useful to investigate suspected cases of poor quality index reads or clashing sequences.

  # -ex/--extension: Allows users to specify the file extension of the demultiplexed FASTQ files produced.

# To see all options, run: python Demultiplexor.py -h


##################
##### OUTPUT #####  
##################
    
# A fastq file will be produced for each sample listed in the index file, in the modified format, containing all reads that matched
# So we go from:        R1 - [6s|X1|----J(D)V-----]  
#                       R2 - [X2]
#                       R3 - [8s|N1|8s|N2|2s|-----5'UTR-----]
# To: ========>         out- [8s|N1|8s|N2|2s|X1|X2|----J(D)V-----]         
# Where X = hexamer index sequence, N = random barcode sequence, and Ys = spacer sequences of length Y
  # The 8s sequences can be used downstream to identify the presence and location of random barcode N sequences
  # 2s is also kept to allow for the possibility finding N sequences produced from slipped reads

##################
#### PACKAGES ####  
##################

from __future__ import division
from itertools import izip
import time
import sys
import argparse
import gzip
import os
import Levenshtein as lev
import collections as coll

__version__ = '2.1'

##########################################################
############# READ IN COMMAND LINE ARGUMENTS #############
##########################################################

def args():
  """args(): Obtains command line arguments which dictate the script's behaviour"""

  # Help flag
  parser = argparse.ArgumentParser(
      description='Script to demultiplex FASTQ data produced using the Chain labs\' ligation TCRseq protocol. Please see https://innate2adaptive.github.io/Decombinator/ for details.')
  # Add arguments
  parser.add_argument(
      '-r1', '--read1', type=str, help='Read 1 FASTQ file', required=True)
  parser.add_argument(
      '-r2', '--read2', type=str, help='Read 2 FASTQ file', required=True)
  parser.add_argument(
      '-i1', '--index1', type=str, help='Index read FASTQ file', required=True)
  parser.add_argument(
      '-ix', '--indexlist', type=str, help='File containing sample/index table', required=False)
  parser.add_argument(
      '-t', '--threshold', type=int, help='Edit distance allowed for fuzzy-string index matching (default=2)', required=False, default=2)
  parser.add_argument(
      '-s', '--suppresssummary', action='store_true', help='Output summary data', required=False)
  parser.add_argument(
      '-a', '--outputall', action='store_true', help='Output all possible index combinations', required=False)
  parser.add_argument(
      '-dz', '--dontgzip', action='store_true', help='Stop the output FASTQ files automatically being compressed with gzip (True/False)', required=False)
  parser.add_argument(
      '-dc', '--dontcount', action='store_true', help='Show the count (True/False)', required=False)
  parser.add_argument(
      '-fl', '--fuzzylist', type=bool, help='Output a list of those reads that demultiplexed using fuzzy index matching (True/False)', required=False, default=False)  
  parser.add_argument(
      '-ex', '--extension', type=str, help='Specify the file extension of the output FASTQ files. Default = \"fq\"', required=False, default="fq")
  return parser.parse_args()

############################################
############# FASTQ PROCESSING #############
############################################

def fastq_check(infile):
  """fastq_check(file): Performs a rudimentary sanity check to see whether a file is indeed a FASTQ file"""
  
  success = True
    
  if infile.endswith('.gz'):
    with gzip.open(infile) as possfq:
      read = [next(possfq) for x in range(4)]
  else:
    with open(infile) as possfq:
      read = [next(possfq) for x in range(4)]    
  
  # @ check
  if read[0][0] <> "@":
    success = False
  # Descriptor check
  if read[2][0] <> "+":
    success = False
  # Read/quality match check
  if len(read[1]) <> len(read[3]):
    success = False  
  
  return(success)

def readfq(fp): 
    """readfq(file):Heng Li's Python implementation of his readfq function 
    See https://github.com/lh3/readfq/blob/master/readfq.py"""
    
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def sort_permissions(fl):
  # Need to ensure proper file permissions on output data
    # If users are running pipeline through Docker might otherwise require root access
  if oct(os.stat(fl).st_mode)[4:] != '666':
    os.chmod(str(fl), 0o666)

inputargs = vars(args())

if inputargs['outputall'] == False and not inputargs['indexlist']:
  print "No indexing file provided, and output all option not enabled; one (or both) is required."
  sys.exit()

for f in [inputargs['read1'], inputargs['read2'], inputargs['index1']]:
  if fastq_check(f) == False:
    print "FASTQ sanity check failed reading", f, "- please ensure that this file is properly formatted and/or gzip compressed."
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

outputreads = coll.Counter()
outputreads["Undetermined"] = 0

usedindexes = coll.defaultdict(list)       # This keeps a track of all files that have been generated to house demultiplexed reads

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
fuzzy_count = 0         # number of sequences that were demultiplexed using non-exact index matches
clash_count = 0         # number of fuzzy ID clashes

fuzzies = []            # list to record IDs matched using fuzzy indexes

t0 = time.time() # Begin timer
  
##########################################################
########### LOOP THROUGH ALL READ FILES IN SYNC ##########
######## PROCESS INTO CORRECT FORMAT & DEMULTIPLEX #######
##########################################################

print "Reading input files..."

# Open read files
if inputargs['read1'].endswith('.gz'):
  fq1 = readfq(gzip.open(inputargs['read1']))
else:
  fq1 = readfq(open(inputargs['read1']))
  
if inputargs['index1'].endswith('.gz'):
  fq2 = readfq(gzip.open(inputargs['index1']))
else:
  fq2 = readfq(open(inputargs['index1']))


if inputargs['read2'].endswith('.gz'):
  fq3 = readfq(gzip.open(inputargs['read2']))
else:
  fq3 = readfq(open(inputargs['read2']))

print "Demultiplexing data..."

for record1, record2, record3 in izip(fq1, fq2, fq3):
  # Readfq function with return each read from each file as a 3 part tuple
    # ('ID', 'SEQUENCE', 'QUALITY')
  count += 1  

  if count % 100000 == 0 and inputargs['dontcount'] == False:
    print '\t read', count
  
### NB For non-standard Illumina encoded fastqs, might need to change which fields are carried into fq_* vars
  
  fq_id = record1[0]

  # N relates to barcode random nucleotides, X denotes index bases
  
  ### FORMATTING OUTPUT READ ###
  
  Nseq = record3[1][0:30]
  X1seq = record1[1][6:12]
  X2seq = record2[1]
  readseq = record1[1][12:]

  Nqual = record3[2][0:30]
  X1qual = record1[2][6:12]
  X2qual = record2[2]
  readqual = record1[2][12:]
    
  fq_seq = Nseq + X1seq + X2seq + readseq
  fq_qual = Nqual + X1qual + X2qual + readqual
  
  new_record = str("@" + fq_id + "\n" + fq_seq + "\n+\n" + fq_qual + "\n")  
  
  ### DEMULTIPLEXING ###
  
  seqX = X1seq + X2seq 
  
  if seqX in XXdict:
    # Exact index matches
          
    XXdict[seqX].write(new_record)
    
    dmpd_count += 1
    outputreads[str(XXdict[seqX]).split("\'")[1][:-len(suffix)]] += 1
    
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
      outputreads[str(XXdict[matches[0]]).split("\'")[1][:-len(suffix)]] += 1
      
    else:
      
      if len(matches) > 1:
        clash_count += 1
        
      failed.write(new_record)
      outputreads['Undetermined'] += 1
  
for x in XXdict.values():
  x.close()
  sort_permissions(x.name)

failed.close()
sort_permissions(failed.name)
fq1.close()
fq2.close()
fq3.close()

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
        sort_permissions(outfile.name)
    os.unlink(f + suffix)

#################################################
################## STATISTICS ###################
#################################################

timed = time.time() - t0
took = round(timed,2)
#print count, 'reads processed from', rd1file, 'and', fq2file, 'and output into', outfq #FIX
if took < 60:
  print '\t\t\t\t\t\t\tTook', took, 'seconds to demultiplex samples'
else:
  print '\t\t\t\t\t\t\tTook', round((timed/60),2), 'minutes to jimmy indexes and hexamers around'


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
  
  # Check for existing date-stamped file
  summaryname = "Logs/" + date + "_Demultiplexing_Summary.csv"
  if not os.path.exists(summaryname): 
    summaryfile = open(summaryname, "w")
  else:
    # If one exists, start an incremental day stamp
    for i in range(2,10000):
      summaryname = "Logs/" + date + "_Demultiplexing_Summary_" + str(i) + ".csv"
      if not os.path.exists(summaryname): 
        summaryfile = open(summaryname, "w")
        break
      
  # Generate string to write to summary file
  
  summstr = "Property,Value\nDirectory," + os.getcwd() + "\nDateFinished," + date + "\nTimeFinished," + time.strftime("%H:%M:%S") + "\nTimeTaken(Seconds)," + str(round(timed,2)) + "\n"
  
  for s in ['read1', 'read2', 'index1', 'indexlist', 'extension', 'threshold', 'outputall', 'dontgzip', 'fuzzylist']:
    summstr = summstr + s + "," + str(inputargs[s]) + "\n"
  
  summstr = summstr + "NumberReadsInput," + str(count) + "\nNumberReadsDemultiplexed," + str(dmpd_count) + "\nNumberFuzzyDemultiplexed," + str(fuzzy_count) + "\nNumberIndexClash," + str(clash_count) + "\n\nOutputFile,IndexUsed\n" 
  
    # Write out number of reads in and details of each individual output file
  for x in sorted(usedindexes.keys()):
    summstr = summstr + x + "," + usedindexes[x] + "\n"
  
  if inputargs['indexlist']:
    summstr = summstr + "\nOutputFile,IndexNumbersUsed(SP1&SP2)\n"
    for x in indexes:
      splt = x.rstrip().split(",")
      summstr = summstr + splt[0] + "," + splt[1] + " & " + splt[2] + "\n"
  
  
  summstr = summstr + "\nOutputFile,NumberReads\n"
  
  for x in sorted(outputreads.keys()):
    summstr = summstr + x + "," + str(outputreads[x]) + "\n"
    
  print >> summaryfile, summstr 
  
  summaryfile.close()
  sort_permissions(summaryname)
  
# Write out list of fuzzy matched sequences, so can fish out later if needed
if inputargs['threshold'] > 0 and inputargs['fuzzylist'] == True:
  
  print "\nOutputting list of reads demultiplexed using fuzzy index matching"
  
  # Check for directory and make summary file
  if not os.path.exists('Logs'):
    os.makedirs('Logs')
  date = time.strftime("%Y_%m_%d")
  
  fuzzname = "Logs/" + date + "_FuzzyMatchedIDs.txt"
  fuzzout = open(fuzzname, "w")
  for f in fuzzies:
    print >> fuzzout, f

  fuzzout.close()
  sort_permissions(fuzzname)
  

