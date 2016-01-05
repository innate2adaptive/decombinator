# FIX - remove "print v/j" etc, and then find where the missing unassigned are going
import traceback    
import sys          
import os
import urllib2
import numpy as np
import decimal as dec
import string
import operator as op
import collections as col
import argparse
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from acora import AcoraBuilder
from time import time, strftime
from string import Template
from operator import itemgetter, attrgetter
import Levenshtein as lev

# ADD
  # species/gamma delta

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
      '-fq', '--fastq', type=str, help='Correctly demultiplexed/processed FASTQ file containing TCR reads', required=True)
  parser.add_argument(
      '-c', '--chain', type=str, help='TCR chain (a/b/g/d)', required=False)
  parser.add_argument(
      '-s', '--suppresssummary', type=bool, help='Output summary data (True/False)', required=False, default=False)
  parser.add_argument(
      '-dz', '--dontgzip', type=bool, help='Stop the output FASTQ files automatically being compressed with gzip (True/False)', required=False, default=False)
  parser.add_argument(
      '-dc', '--dontcount', type=bool, help='Show the count (True/False)', required=False, default=False)
  parser.add_argument(
      '-ex', '--extension', type=str, help='Specify the file extension of the output DCR file. Default = \"n12\"', required=False, default="n12")
  parser.add_argument(
      '-fr', '--frames', type=str, help='Specify the frames to search in (forward/reverse/both). Default = reverse', required=False, default="reverse")
  parser.add_argument(
      '-tg', '--tags', type=str, help='Specify which Decombinator tag set to use (extended or original). Default = extended', required=False, default="extended")
  parser.add_argument(
      '-sp', '--species', type=str, help='Specify which species TCR repertoire the data consists of (human or mouse). Default = human', required=False, default="human")
  parser.add_argument(
      '-N', '--allowNs', type=bool, help='Whether to allow VJ rearrangements containing ambiguous base calls (\'N\'). Default = False', required=False, default=False)
  parser.add_argument(
      '-ln', '--lenthreshold', type=int, help='Acceptable threshold for inter-tag (V to J) sequence length. Default = 130', required=False, default=130)

  return parser.parse_args()


##########################################################
############# FASTQ SANITY CHECK AND PARSING #############
##########################################################

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

def revcomp(read):
  """rc(read): Wrapper for SeqIO reverse complement function"""
  return str(Seq(read).reverse_complement())


def readfq(fp): # this is a generator function
    """readfq(file):Heng Li's Python implementation of his readfq function 
    https://github.com/lh3/readfq/blob/master/readfq.py"""
    
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



#####################################
############# DECOMBINE #############
#####################################

def vanalysis(read):
  
  hold_v = v_key.findall(read)

  if hold_v:
    if len(hold_v) > 1:
      counts['multiple_v_matches'] += 1
      return

    v_match = v_seqs.index(hold_v[0][0]) # Assigns V
    temp_end_v = hold_v[0][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
    
    v_seq_start = hold_v[0][1]      
    end_v_v_dels = get_v_deletions( read, v_match, temp_end_v, v_regions )      
    if end_v_v_dels: # If the number of deletions has been found
      print "v"
      return v_match, end_v_v_dels[0], end_v_v_dels[1], v_seq_start
      
  else:
    
    hold_v1 = half1_v_key.findall(read)
    
    if hold_v1:
      for i in range(len(hold_v1)):
        indices = [y for y, x in enumerate(half1_v_seqs) if x == hold_v1[i][0] ]
        for k in indices:
          if len(v_seqs[k]) == len(read[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[half1_v_seqs.index(hold_v1[i][0])])]):
            if lev.hamming( v_seqs[k], read[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[k])] ) <= 1:
              counts['verr2'] += 1
              v_match = k
              temp_end_v = hold_v1[i][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
              end_v_v_dels = get_v_deletions( read, v_match, temp_end_v, v_regions )
              v_seq_start = hold_v1[i][1]  
              print "v"
              return v_match, end_v_v_dels[0], end_v_v_dels[1], v_seq_start
    
    else:
      hold_v2 = half2_v_key.findall(read)
      if hold_v2:
        for i in range(len(hold_v2)):
          indices = [y for y, x in enumerate(half2_v_seqs) if x == hold_v2[i][0] ]
          for k in indices:
            if len(v_seqs[k]) == len(read[hold_v2[i][1]-v_half_split:hold_v2[i][1]-v_half_split+len(v_seqs[half2_v_seqs.index(hold_v2[i][0])])]):
              if lev.hamming( v_seqs[k], read[hold_v2[i][1]-v_half_split:hold_v2[i][1]+len(v_seqs[k])-v_half_split] ) <= 1:
                counts['verr1'] += 1
                v_match = k
                temp_end_v = hold_v2[i][1] + jump_to_end_v[v_match] - v_half_split - 1 # Finds where the end of a full V would be
                end_v_v_dels = get_v_deletions( read, v_match, temp_end_v, v_regions )
                v_seq_start = hold_v2[i][1] - v_half_split      
                print "v"
                return v_match, end_v_v_dels[0], end_v_v_dels[1], v_seq_start
      else:
        counts['no_v_assigned'] += 1
        return
      
def janalysis(read):
  
  hold_j = j_key.findall(read)
  
  if hold_j:
    if len(hold_j) > 1:
      counts['multiple_j_matches'] += 1
      return
  
    j_match = j_seqs.index(hold_j[0][0]) # Assigns J
    temp_start_j = hold_j[0][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
    
    j_seq_end = hold_j[0][1] + len(hold_j[0][0])      
        
    start_j_j_dels = get_j_deletions( read, j_match, temp_start_j, j_regions )
    
    if start_j_j_dels: # If the number of deletions has been found
      print "j"
      return j_match, start_j_j_dels[0], start_j_j_dels[1], j_seq_end
          
  else:
    
    hold_j1 = half1_j_key.findall(read)
    if hold_j1:
      for i in range(len(hold_j1)):
        indices = [y for y, x in enumerate(half1_j_seqs) if x == hold_j1[i][0] ]
        for k in indices:
          if len(j_seqs[k]) == len(read[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[half1_j_seqs.index(hold_j1[i][0])])]):
            if lev.hamming( j_seqs[k], read[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[k])] ) <= 1:
              counts['jerr2'] += 1
              j_match = k
              temp_start_j = hold_j1[i][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
              j_seq_end = hold_j1[i][1] + len(hold_j1[i][0]) + j_half_split                                              
              start_j_j_dels = get_j_deletions( read, j_match, temp_start_j, j_regions )
              print "j"
              return j_match, start_j_j_dels[0], start_j_j_dels[1], j_seq_end
            
    else:        
      hold_j2 = half2_j_key.findall(read)
      if hold_j2:
        for i in range(len(hold_j2)):
          indices = [y for y, x in enumerate(half2_j_seqs) if x == hold_j2[i][0] ]
          for k in indices:
            if len(j_seqs[k]) == len(read[hold_j2[i][1]-j_half_split:hold_j2[i][1]-j_half_split+len(j_seqs[half2_j_seqs.index(hold_j2[i][0])])]):
              if lev.hamming( j_seqs[k], read[hold_j2[i][1]-j_half_split:hold_j2[i][1]+len(j_seqs[k])-j_half_split] ) <= 1:
                counts['jerr1'] += 1
                j_match = k
                temp_start_j = hold_j2[i][1] - jump_to_start_j[j_match] - j_half_split # Finds where the start of a full J would be
                j_seq_end = hold_j2[i][1] + len(hold_j2[i][0])                                                
                start_j_j_dels = get_j_deletions( read, j_match, temp_start_j, j_regions )
                print "j"
                return j_match, start_j_j_dels[0], start_j_j_dels[1], j_seq_end
      else:
         counts['no_j_assigned'] += 1
         return
       
def dcr(read):

  """dcr(read): Core function which checks a read (in the given frame) for a rearranged TCR of the specified chain.
    Returns a list giving: V gene index, J gene index, # deletions in V gene, # deletions in J gene,
      insert sequence (between ends of V and J), inter-tag sequence (for collapsing), and its quality scores"""

  v_seq_start = 0     
  j_seq_end = 0      
  
  vdat = vanalysis(read)
  
  if not vdat:
    #counts['VTEST'] += 1 #FIXDEL
    return

  jdat = janalysis(read)
  
  if jdat:
    
    # Filter out rearrangements with indications they probably represent erroneous sequences
    if "N" in read[vdat[3]:jdat[3]] and inputargs['allowNs'] == False:                          # Ambiguous base in inter-tag region
      counts['dcrfilter_intertagN'] += 1
    elif (vdat[3] - jdat[3]) >= inputargs['lenthreshold']:                                      # Inter-tag length threshold
      counts['dcrfilter_toolong_intertag'] += 1
    elif vdat[2] > (jump_to_end_v[vdat[0]] - len(v_seqs[vdat[0]])) or jdat[2] > jump_to_start_j[jdat[0]]: # Impossible number of deletions
      counts['dcrfilter_imposs_deletion'] += 1                    
    elif (vdat[3] + len(v_seqs[vdat[0]])) > (jdat[3] + len(j_seqs[jdat[0]])):                             # Overlapping tags 
      counts['dcrfilter_tag_overlap'] += 1                     
    
    else:        
      vj_details = [vdat[0], jdat[0], vdat[2], jdat[2], read[vdat[1]+1:jdat[1]], vdat[3], jdat[3]]
      return vj_details
  
  else:
    counts['VJ_assignment_failed'] += 1
    return
  

###########################################################
############# ANCILLARY DECOMBINING FUNCTIONS #############
###########################################################


def get_v_deletions( read, v_match, temp_end_v, v_regions_cut ):
    # This function determines the number of V deletions in sequence read
    # by comparing it to v_match, beginning by making comparisons at the
    # end of v_match and at position temp_end_v in read.
    function_temp_end_v = temp_end_v
    pos = -1
    is_v_match = 0
    
    # Catch situations in which the temporary end of the V exists beyond the end of the read
    if function_temp_end_v >= len(read):
      counts['v_del_failed_tag_at_end'] += 1
      return
      
    while is_v_match == 0 and 0 <= function_temp_end_v < len(read):
        if str(v_regions_cut[v_match])[pos] == read[function_temp_end_v] and str(v_regions_cut[v_match])[pos-1] == read[function_temp_end_v-1] and str(v_regions_cut[v_match])[pos-2] == read[function_temp_end_v-2]:
            is_v_match = 1
            deletions_v = -pos - 1
            end_v = function_temp_end_v
        else:
            pos -= 1
            function_temp_end_v -= 1

    if is_v_match == 1:
        return [end_v, deletions_v]
    else:
        counts['v_del_failed'] += 1
        return 

def get_j_deletions( read, j_match, temp_start_j, j_regions_cut ):
    # This function determines the number of J deletions in sequence read
    # by comparing it to j_match, beginning by making comparisons at the
    # end of j_match and at position temp_end_j in read.
    function_temp_start_j = temp_start_j
    pos = 0
    is_j_match = 0
    
    while is_j_match == 0 and 0 <= function_temp_start_j+2 < len(str(read)):
        if str(j_regions_cut[j_match])[pos] == read[function_temp_start_j] and str(j_regions_cut[j_match])[pos+1] == read[function_temp_start_j+1] and str(j_regions_cut[j_match])[pos+2] == read[function_temp_start_j+2]:
            is_j_match = 1
            deletions_j = pos
            start_j = function_temp_start_j
        else:
            pos += 1
            function_temp_start_j += 1
            
    if is_j_match == 1:
        return [start_j, deletions_j]
    else:
        counts['j_del_failed'] += 1
        return 

def get_v_tags(file_v, half_split):
    import string
    
    v_seqs = [] # Holds all V tags
    jump_to_end_v = [] # Holds the number of jumps to make to look for deletions for each V region once the corresponding tag has been found
    for line in file_v:
        elements = line.rstrip("\n") # Turns every element in a text file line separated by a space into elements in a list
        v_seqs.append(string.split(elements)[0]) # Adds elements in first column iteratively
        jump_to_end_v.append(int(string.split(elements)[1])) # Adds elements in second column iteratively

    half1_v_seqs = []
    half2_v_seqs = []

    for i in range(len(v_seqs)):
        half1_v_seqs.append(v_seqs[i][0:half_split])
        half2_v_seqs.append(v_seqs[i][half_split:])
    
    return [v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v]

def get_j_tags(file_j, half_split):
    import string
    
    j_seqs = [] # Holds all J tags
    jump_to_start_j = [] # Holds the number of jumps to make to look for deletions for each J region once the corresponding tag has been found

    for line in file_j:
        elements = line.rstrip("\n")
        j_seqs.append(string.split(elements)[0])
        jump_to_start_j.append(int(string.split(elements)[1]))

    half1_j_seqs = []
    half2_j_seqs = []

    for j in range(len(j_seqs)):
        half1_j_seqs.append(j_seqs[j][0:half_split])
        half2_j_seqs.append(j_seqs[j][half_split:])

    return [j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j]



##########################################################
############# READ IN COMMAND LINE ARGUMENTS #############
##########################################################


inputargs = vars(args())

counts = col.Counter()
counts['start_time'] = time()

# Brief FASTQ sanity check
if fastq_check(inputargs['fastq']) <> True:
  print "FASTQ sanity check failed reading", inputargs['fastq'], "- please ensure that this file is properly formatted and/or gzip compressed."
  sys.exit()

# Get chain information
if inputargs['chain'].upper() in ['A', 'ALPHA', 'TRA', 'TCRA']:
  chain = "a" 
elif inputargs['chain'].upper() in ['B', 'BETA', 'TRB', 'TCRB']:
  chain = "b" 
elif inputargs['chain'].upper() in ['G', 'GAMMA', 'TRG', 'TCRG']:
  chain = "g" 
elif inputargs['chain'].upper() in ['D', 'DELTA', 'TRD', 'TCRD']:
  chain = "d" 
else:
  print "TCR chain not recognised. Please choose from a/b/g/d (case-insensitive)."
  sys.exit()


#################################################
############# GET GENES, BUILD TRIE #############
#################################################

chainnams = {"a": "alpha", "b": "beta", "g": "gamma", "d": "delta"}

print 'Importing', chainnams[chain], 'TCR gene sequences...'

# Do not change - V tags are split at position 10, J at position 10, to look for half tags if no full tag is found.
v_half_split, j_half_split = [10,10] 

# FIX - need to use appropriate tag splitting options depending on tag set (which is also species dependent)

  
# Look for tag and V/J fasta and tag files: if these cannot be found in the working directory, source them from GitHub repositories
if inputargs['tags'] == "extended" and inputargs['species'] == "mouse":
  print "There is currently not an extended tag set for mouse TCR genes.\n \
  Please use the following flags: -sp mouse -tg original"
  sys.exit()

if inputargs['tags'] == "extended":
  for gene in ['v', 'j']:
    local_fasta = "exthuman_TR" + chain.upper() + gene.upper() + "_region.fasta"
    if os.path.isfile(local_fasta):
      with open(local_fasta) as infile:
        vars()[gene + "_genes"] = list(SeqIO.parse(infile, "fasta"))
    else:
      online_fasta = urllib2.urlopen("https://raw.githubusercontent.com/JamieHeather/tcr-analysis/master/exthuman_TR" + chain.upper() + gene.upper() + "_region.fasta")
      vars()[gene + "_genes"] = list(SeqIO.parse(online_fasta, "fasta"))
      online_fasta.close()
        
elif inputargs['tags'] == "original":  
  for gene in ['v', 'j']:
    local_fasta = inputargs['species'] + "_TR" + chain.upper() + gene.upper() + "_region.fasta"
    if os.path.isfile(local_fasta):
      with open(local_fasta) as infile:
        vars()[gene + "_genes"] = list(SeqIO.parse(iinfile, "fasta"))
    else:
      online_fasta = urllib2.urlopen("https://raw.githubusercontent.com/uclinfectionimmunity/Decombinator/master/" + inputargs['species'] \
        + "_TR" + chain.upper() + gene.upper() + "_region.fasta")
      vars()[gene + "_genes"] = list(SeqIO.parse(online_fasta, "fasta"))
      online_fasta.close()
      
else:
  print "Tag set unrecognised. Check tag set and species flag."
  sys.exit()

v_regions = []
for j in range(0, len(v_genes)):
    v_regions.append(string.upper(v_genes[j].seq))

j_regions = []
for j in range(0, len(j_genes)):
    j_regions.append(string.upper(j_genes[j].seq))



##############
## Build keyword tries of V and J inputargs['tags'] for fast assignment
  # FIX need to sort this out to be species and gene specific
if inputargs['tags'] == "original":
  if os.path.isfile("tags_tr" + chain.lower() + "v.txt") and os.path.isfile("tags_tr" + chain.lower() + "j.txt"):
    v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(open("tags_tr" + chain.lower() + "v.txt", "rU"), v_half_split)
    j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(open("tags_tr" + chain.lower() + "j.txt", "rU"), j_half_split)
  else:
    v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(urllib2.urlopen("https://raw.githubusercontent.com/uclinfectionimmunity/Decombinator/master/humantags_tr" + chain.lower() + "v.txt"), v_half_split)
    j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(urllib2.urlopen("https://raw.githubusercontent.com/uclinfectionimmunity/Decombinator/master/humantags_tr" + chain.lower() + "j.txt"), j_half_split)               
elif inputargs['tags'] == "extended":
  if os.path.isfile("exttags_tr" + chain.lower() + "v.txt") and os.path.isfile("exttags_tr" + chain.lower() + "j.txt"):  
    v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(open("exttags_tr" + chain.lower() + "v.txt", "rU"), v_half_split)
    j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(open("exttags_tr" + chain.lower() + "j.txt", "rU"), j_half_split)
  else:
    v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(urllib2.urlopen("https://raw.githubusercontent.com/JamieHeather/tcr-analysis/master/exttags_tr" + chain.lower() + "v.txt"), v_half_split)
    j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(urllib2.urlopen("https://raw.githubusercontent.com/JamieHeather/tcr-analysis/master/exttags_tr" + chain.lower() + "j.txt"), j_half_split) 

v_builder = AcoraBuilder()
for i in range(0,len(v_seqs)):
    v_builder.add(str(v_seqs[i])) # Add all V tags to keyword trie

v_key = v_builder.build()

j_builder = AcoraBuilder()
for i in range(0,len(j_seqs)):
    j_builder.add(str(j_seqs[i])) # Add all J tags to keyword trie

j_key = j_builder.build()

##############
## Build keyword tries for first and second halves of both V and J tags
v_half1_builder = AcoraBuilder()
for i in range(0,len(half1_v_seqs)):
    v_half1_builder.add(str(half1_v_seqs[i]))
half1_v_key = v_half1_builder.build()

v_half2_builder = AcoraBuilder()
for i in range(0,len(half2_v_seqs)):
    v_half2_builder.add(str(half2_v_seqs[i]))
half2_v_key = v_half2_builder.build()

j_half1_builder = AcoraBuilder()
for i in range(0,len(half1_j_seqs)):
    j_half1_builder.add(str(half1_j_seqs[i]))
half1_j_key = j_half1_builder.build()

j_half2_builder = AcoraBuilder()
for i in range(0,len(half2_j_seqs)):
    j_half2_builder.add(str(half2_j_seqs[i]))
half2_j_key = j_half2_builder.build()



#########################################################
############# SCROLL THROUGH FILE & ANALYSE #############
#########################################################

print "Decombining FASTQ data..."

suffix = "." + inputargs['extension']
samplenam = str(inputargs['fastq'].split(".")[0]) 

if inputargs['tags'] == "original":
  name_results = "vDCR_" + chainnams[chain] + "_" + samplenam
elif inputargs['tags'] == "extended":
  name_results = "V3.2_" + chainnams[chain] + "_" + samplenam
else:
  print "Name of tag set not recognised. Please edit code so that -tg = either \'original\' or \'extended\'"
  sys.exit()

stemplate = Template('$v $j $del_v $del_j $nt_insert $seqid $tcr_seq $tcr_qual $barcode $barqual')      

with open(name_results + suffix, 'w') as outfile:

  if inputargs['fastq'].endswith('.gz'):
    opener = gzip.open
  else:
    opener = open
    
  with opener(inputargs['fastq']) as f:
    
    for readid, seq, qual in readfq(f):

      start_time = time()
      vdj = seq[30:]
      
      counts['read_count'] += 1
      if counts['read_count'] % 100000 == 0 and inputargs['dontcount'] == False:
        print '\t read', counts['read_count'] 
  
      # Get details of the VJ recombination
      if inputargs['frames'] == 'reverse':
        recom = dcr(revcomp(vdj))
        frame = 'reverse'
      elif inputargs['frames'] == 'forward':
        recom = dcr(vdj)
        frame = 'forward'
      elif inputargs['frames'] == 'both':
        recom = dcr(revcomp)
        frame = 'reverse'
        if not recom:
          recom = dcr(vdj)
          frame = 'forward'
                      
      if recom:
        counts['vj_count'] += 1
        
        bc = seq[:30]  
        bcQ = qual[:30]
        
        vdjqual = qual[30:]
      
        if frame == 'reverse':
          tcrseq = revcomp(vdj)[recom[5]:recom[6]]
          tcrQ = vdjqual[::-1][recom[5]:recom[6]]
        elif frame == 'forward':
          tcrseq = vdj[recom[5]:recom[6]]
          tcrQ = vdjqual[recom[5]:recom[6]]
        
              
        dcr_string = stemplate.substitute( v = str(recom[0]) + ',', j = str(recom[1]) + ',', del_v = str(recom[2]) + ',', \
          del_j = str(recom[3]) + ',', nt_insert = recom[4] + ',', seqid = readid + ',', tcr_seq = tcrseq + ',', \
          tcr_qual = tcrQ + ',', barcode = bc + ',', barqual = bcQ )      
                              
        outfile.write(dcr_string + '\n')
      
      #times.append(time() - start_time)
      

counts['end_time'] = time()
timetaken = counts['end_time']-counts['start_time']

if inputargs['dontgzip'] == False:
  print "Compressing Decombinator output file,", name_results + suffix, "..."
  
  with open(name_results + suffix) as infile, gzip.open(name_results + suffix + '.gz', 'wb') as outfile:
      outfile.writelines(infile)
  os.unlink(name_results + suffix)

  outfilenam = name_results + suffix + ".gz"
else:
  outfilenam = name_results + suffix


##############################################
############# WRITE SUMMARY DATA #############
##############################################
   

print "Analysed", str(counts['read_count']), "reads, finding", str(counts['vj_count']), chainnams[chain], "VJ rearrangements"
print "Reading from", inputargs['fastq'] + ", writing to", outfilenam
print "Took", str(round(timetaken,2)), "seconds"
#print np.mean(times)



# Write data to summary file
if inputargs['suppresssummary'] == False:
  
  # Check for directory and make summary file
  if not os.path.exists('Logs'):
    os.makedirs('Logs')
  date = strftime("%Y_%m_%d")
  
  # Check for existing date-stamped file
  summaryname = "Logs/" + date + "_" + samplenam + "_Decombinator_Summary.csv"
  if not os.path.exists(summaryname): 
    summaryfile = open(summaryname, "w")
  else:
    # If one exists, start an incremental day stamp
    for i in range(2,10000):
      summaryname = "Logs/" + date + "_" + samplenam + "_Decombinator_Summary" + str(i) + ".csv"
      if not os.path.exists(summaryname): 
        summaryfile = open(summaryname, "w")
        break
      
  # Generate string to write to summary file
  
  summstr = "Property,Value\nDirectory," + os.getcwd() + "\nInputFile," + inputargs['fastq'] + "\nOutputFile," + outfilenam \
    + "\nDateFinished," + date + "\nTimeFinished," + strftime("%H:%M:%S") + "\nTimeTaken(Seconds)," + str(round(timetaken,2)) + "\n\n"
  
  for s in ['species', 'chain','extension', 'tags', 'dontgzip', 'allowNs', 'frames', 'lenthreshold']:
    summstr = summstr + s + "," + str(inputargs[s]) + "\n"
  
  summstr = summstr + "\nNumberReadsInput," + str(counts['read_count']) + "\nNumberReadsDecombined," + str(counts['vj_count']) 
 
  ######################### FIX - MAKE AN EXTENDED OUTPUT OPTION FOR THE FOLLOWING?
  
  # Half tag matching details
  summstr = summstr + "\n\nReadsAssignedUsingHalfTags:,\nV1error," + str(counts['verr1']) \
    + "\nV2error," + str(counts['verr2']) \
    + "\nJ1error," + str(counts['jerr1']) \
    + "\nJ2error," + str(counts['jerr2'])

  
  # Number reads filtered out
  summstr = summstr + "\n\nReadsFilteredOut:,\nAmbiguousBaseCall," + str(counts['dcrfilter_intertagN']) \
    + "\nOverlongInterTagSeq," + str(counts['dcrfilter_toolong_intertag']) \
    + "\nImpossibleDeletions," + str(counts['dcrfilter_imposs_deletion']) \
    + "\nOverlappingTagBoundaries," + str(counts['dcrfilter_tag_overlap']) \
      
  ##########################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##########################
  # NEEDS FIXING - sum of all failures does not match the number of reads that failed to decombine, need to figure out why
  # Number reads failed to assign 
  #summstr = summstr + "\n\nReadsFailedAssignment:,\nMultipleVtagMatches," + str(counts['multiple_v_matches']) \
    #+ "\nVTagAtEndRead," + str(counts['v_del_failed_tag_at_end']) \
    #+ "\nVDeletionsUndetermined," + str(counts['v_del_failed']) \
    #+ "\nNoVDetected," + str(counts['no_v_assigned']) \
    #+ "\nMultipleJTagMatches," + str(counts['multiple_j_matches']) \
    #+ "\nJDeletionsUndermined," + str(counts['j_del_failed']) \
    #+ "\nNoJDetected," + str(counts['no_j_assigned']) \
    #+ "\nVJGeneAssignmentFailed," + str(counts['VJ_assignment_failed'])     
      
  print >> summaryfile, summstr 
  
  summaryfile.close()
  




sys.exit()



