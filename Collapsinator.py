# Collapsinator
# Katharine Best and James M. Heather, August 2016, UCL
# https://innate2adaptive.github.io/Decombinator/

##################
### BACKGROUND ###
##################

# Takes the output files of Decombinator (run using the barcoding option) and performs collapsing and error correction
# This version is a modified version of KB's script collapsinator_20141126.py
  # That was itself an improved version of the CollapseTCRs.py script used in the Heather et al HIV TCR paper (DOI: 10.3389/fimmu.2015.00644)
# Version 4.0.1 includes an improved statistical clustering method by Peter de Greef, Utrecht University

##################
###### INPUT #####
##################

# Takes as input .n12 files produced by Decombinator (v3), assuming it has been run on suitably barcoded and demultiplexed data.

# Other optional flags:
  
  # -s/--supresssummary: Supress the production of a summary file containing details of the run into a 'Logs' directory. 
  
  # -dz/--dontgzip: Suppress the automatic compression of output demultiplexed FASTQ files with gzip. 
  
  # -dc/--dontcount: Suppress the whether or not to show the running line count, every 100,000 reads. 
    # Helps in monitoring progress of large batches. D

  # The other optional flags are somewhat complex, and caution is advised in their alteration.

# To see all options, run: python Collapsinator.py -h

# Input files need to be in the appropriate format, consisting of:
  # V index, J index, # V deletions, # J deletions, insert, ID, inter-tag TCR sequence, inter-tag quality, barcode sequence, barcode quality

##################
##### OUTPUT #####  
##################
   
# A Decombinator index file, giving each error-corrected DCR index and the frequency with which it appears in the final processed data

########################################################################################################################
from __future__ import division
import collections as coll
import random
import Levenshtein as lev
from time import time, strftime
import argparse
import gzip
import regex
from scipy.special import comb
import copy
import os, sys

from IPython import embed

__version__ = '4.0.1'

##########################################################
############# READ IN COMMAND LINE ARGUMENTS #############
##########################################################

def args():
  """args(): Obtains command line arguments which dictate the script's behaviour"""

  # Help flag
  parser = argparse.ArgumentParser(
      description='Collapse and error correct Decombined TCR sequences produced using the Chain lab\'s ligation TCRseq protocol. Please see https://innate2adaptive.github.io/Decombinator/ for details.')
  # Add arguments
  parser.add_argument(
      '-in', '--infile', type=str, help='File containing raw verbose Decombinator output, i.e. \
      5 part classifier plus barcode and inter-tag sequence and quality strings', required=True)
  parser.add_argument(
      '-mq', '--minbcQ', type=int, help='Minimum quality score that barcode nucleotides should be to for that rearrangement to be retained. Default = 20.', \
        required=False, default=20)
  parser.add_argument(
      '-bm', '--bcQbelowmin', type=int, help='Number of nucleotides per barcode whose quality score are allowed to be below -mq and still be retained. Default = 1.', \
        required=False, default=1)
  parser.add_argument(
      '-aq', '--avgQthreshold', type=int, help='Average quality threshold that barcode sequences must remain above for rearrangements to be retained. Default = 30', \
        required=False, default=30)
  parser.add_argument(
      '-lv', '--percentlevdist', type=int, help='Percentage Levenshtein distance that is allowed to estimate whether two sequences within a barcode are \
      derived from the same originator molecule. Default = 10', required=False, default=10)  
  parser.add_argument(
      '-bc', '--bcthreshold', type=int, help='Number of sequence edits that are allowed to consider two barcodes to be derived from same originator \
      during clustering. Default = 2.', required=False, default=2)
  parser.add_argument(
      '-s', '--suppresssummary',  action='store_true', help='Suppress the production of output summary data log', required=False)
  parser.add_argument(
      '-dz', '--dontgzip', action='store_true', help='Stop the output FASTQ files automatically being compressed with gzip', required=False)
  parser.add_argument(
      '-dc', '--dontcount', action='store_true', help='Block the line count from being shown', required=False)
  parser.add_argument(
      '-ex', '--extension', type=str, help='Specify the file extension of the output DCR file. Default = \"freq\"', required=False, default="freq")
  parser.add_argument(
      '-N', '--allowNs', action='store_true', help='Used to allow VJ rearrangements containing ambiguous base calls (\'N\')', required=False)
  parser.add_argument(
      '-ln', '--lenthreshold', type=int, help='Acceptable threshold for inter-tag (V to J) sequence length. Default = 130', required=False, default=130)
  parser.add_argument(
      '-di', '--dontcheckinput', action='store_true', help='Override the inputfile sanity check', required=False)
  parser.add_argument(
      '-bd', '--barcodeduplication', action='store_true', help='Optionally output a file containing the final list of clustered barcodes, and their frequencies',\
        required=False)
  parser.add_argument(
      '-pb', '--positionalbarcodes', action='store_true', help='Instead of inferring random barcode sequences from their context relative to spacer sequences, just take the sequence at the default positions. Useful to salvage runs when R2 quality is terrible.',\
        required=False)
  parser.add_argument(
      '-ol', '--oligo', type=str, help='Choose experimental oligo for correct identification of spacers ["M13", "I8"] (default: M13)',\
        required=False, default="m13")
  
  return parser.parse_args()

    
########################################################################################################################
# Functions

def num_check(poss_int):
    """ Check whether string is feasibly an integer value of zero or greater """
    try:
      if int(poss_int) >= 0: 
        return True
      else:
        return False
    except ValueError:
      return False
    
def is_dna(poss_dna):
    """ Check whether string is feasibly a DNA sequence read"""
    return set(poss_dna.upper()).issubset({'A', 'C', 'G', 'T', 'N'})

def flatten(l):
  return [item for sublist in l for item in sublist]
      
def check_dcr_file(infile):
    """ 
    Perform sanity check on input Decombinator file
    Check first few lines to see whether they fit the correct criteria this script relies on
    """
    # Check whether file a) exists and b) is not empty
    if os.path.isfile(infile) == False:
      print('Cannot find file, please double-check path.')
      return False
    if os.path.getsize(infile) == 0:
      print('Input file appears to be empty; please double-check path.')
      return False
    
    # Check first few lines    
    with opener(infile,"rt") as poss_dcr: 
      for i in range(5):
              
        # Check it's a comma-delimited file
        if "," not in next(poss_dcr):
          print('Input Decombinator file sanity check fail: seemingly not comma-delimited text file.')
          return False
        else:
          fields = next(poss_dcr).rstrip().split(", ")
        
        # Check lines contain correct number of fields
        if len(fields) != 10:
          print('Input Decombinator file sanity check fail: file does not contain the correct number of comma-delimited fields (ten).')
          return False
        # Check DCR classifiers are feasible (i.e. 4x integers followed by a DNA string)
        elif not all([num_check(field) for field in fields[0:4]]):
          print('Input Decombinator file sanity check fail: integer components of Decombinator classifier not feasible (i.e. not integers >= zero).')
          return False        
        elif not is_dna(fields[4]):
          print('Input Decombinator file sanity check fail: Decombinator insert sequence field contains non-DNA sequence.')
          return False
        # Check inter-tag and barcode sequences are legitimate DNA strings 
        elif is_dna(fields[6]) != True or is_dna(fields[8]) != True:
          print('Input Decombinator file sanity check fail: inter-tag and/or barcode sequence contains non-DNA sequences.')
          return False
        # Check inter-tag and barcode quality are feasible FASTQ quality scores
        elif all(set([num_check(x) for x in get_qual_scores(fields[7])])) != True or all(set([num_check(x) for x in get_qual_scores(fields[9])])) != True:
          print('Input Decombinator file sanity check fail: inter-tag and/or barcode quality strings do not appear to contain valid scores.')
          return False        
        # Check that inter-tag and barcode sequence/quality pairs are the correct length
        elif len(fields[6]) != len(fields[7]) or len(fields[8]) != len(fields[9]):
          print('Input Decombinator file sanity check fail: inter-tag and/or barcode sequence and quality string pairs are not of the same length.')
          return False    
    
      return True       # If first few lines all pass, assume the file is fine and carry on with the analysis.

def getOligo(oligo_name):
  # New oligos can be added here, specifying their spacers in the given format, and adding them to
  # the returned list.
  oligos = {}
  oligos["m13"] = {"spcr1": "GTCGTGACTGGGAAAACCCTGG","spcr2":"GTCGTGAT"}
  oligos["i8"] = {"spcr1":"GTCGTGAT","spcr2":"GTCGTGAT"}

  if oligo_name.lower() not in oligos:
    print("Error: Failed to recognise oligo name. Please choose from " + str(list(oligos.keys())))
    sys.exit()  
  
  return oligos[oligo_name.lower()]

def findSubs(subseq, seq):
    # allow up to two substitutions in subseq.
    err_subseqs = regex.findall("("+subseq+"){1s<=2}", seq)
    return err_subseqs

def findSubsInsOrDels(subseq, seq):
    # allow up to two subsitutions OR (one deletion or one insertion) in subseq.
    err_subseqs = regex.findall("("+subseq+"){2i+2d+1s<=2}", seq)
    return err_subseqs

def spacerSearch(subseq,seq):
    # first search for subseq within seq. If unsuccessful, allow for subsitutions in subseq.
    # If still unsuccessful, allow for deletions, insertions and subsitutions in subseq.
    foundseq = regex.findall(subseq, seq)
    if not foundseq:
      foundseq = findSubs(subseq, seq)
    if not foundseq:
      foundseq = findSubsInsOrDels(subseq, seq)
    return foundseq

def findFirstSpacer(oligo,seq):
    allowance = 4
    spacer = []
    spcr1 = oligo["spcr1"]
    spacer += spacerSearch(spcr1, seq[0:len(spcr1)+allowance])
    return spacer

def findSecondSpacer(oligo,seq):
    spacer = []
    spcr1 = oligo["spcr1"]
    spcr2 = oligo["spcr2"]
    spacer += spacerSearch(spcr2,seq[len(spcr1):])
    return spacer

def getSpacerPositions(bcseq,spacers):
    # locates the positions of the determined spacers in the sequence
    positions = []
    startpos = 0
    for x in spacers:
      positions.append(bcseq.find(x,startpos))
      startpos += len(x)
    return positions

def filterShortandLongBarcodes(b1len,b2end,bcseq,counts):
    if b1len <= 3:
      counts['getbarcode_fail_n1tooshort'] += 1
      return 'fail'
    elif b1len >= 9:
      counts['getbarcode_fail_n1toolong'] += 1
      return 'fail'
    elif b2end > len(bcseq):
      counts['getbarcode_fail_n2pastend'] += 1
      return 'fail'
    else:
      return 'pass'

def logExactOrRegexMatch(spacers,oligo,counts):
    if spacers == [oligo['spcr1'],oligo['spcr2']]:
      counts['getbarcode_pass_exactmatch'] += 1
    else:
      counts['getbarcode_pass_regexmatch'] += 1

def logFuzzyMatching(b1len,bclength,spacers,oligo,counts):
    if b1len == bclength and spacers != [oligo['spcr1'],oligo['spcr2']]:
      counts['getbarcode_pass_fuzzymatch_rightlen'] += 1  
    elif b1len in [4,5] and spacers != [oligo['spcr1'],oligo['spcr2']]:
      counts['getbarcode_pass_fuzzymatch_short'] += 1 
    elif b1len >= 7 and spacers != [oligo['spcr1'],oligo['spcr2']]:
      counts['getbarcode_pass_fuzzymatch_long'] += 1
    elif b1len == bclength:
       counts['getbarcode_pass_other'] += 1

def get_barcode(bcseq,inputargs,counts):
  """
  Given a barcode-region sequence, outputs the sequence of the do-docamer barcode.
  This barcode (theoretically) consists of the concatentation of the two random hexamer sequences contained in the ligation oligo.
  However errors in sequences and ligation oligo production can mean that the random nucleotides are not always at the expected position.
  This function uses the known sequence of the spacers (which bound each of the two N6s to their 5') to deduce the random sequences.
  Returns a list of four numbers, giving the start and stop positions of N1 and N2 respectively.
  """ 
  if "N" in bcseq and inputargs['allowNs'] == False:    # ambiguous base-call check 
    counts['getbarcode_fail_N'] += 1
    return

  # gets spacer sequences of specified oligo
  oligo = getOligo(inputargs['oligo'])
  
  # sets first spacer based on specified oligo
  spacers = findFirstSpacer(oligo, bcseq)  

  # sequences with no first spacer are removed from analysis
  if not len(spacers) == 1:
    return None

  # sets second spacer based on specified oligo
  spacers += findSecondSpacer(oligo, bcseq)

  # sequences which do not have two spacers are logged then removed from analysis
  if not len(spacers) == 2:
    counts['getbarcode_fail_not2spacersfound'] += 1
    return None

  spacer_positions = getSpacerPositions(bcseq, spacers)

  # set expected barcode length
  bclength = 6
  # start and end of barcode positions are set
  b1start = spacer_positions[0] + len(spacers[0])
  b1end   = spacer_positions[1]
  b2start = spacer_positions[1] + len(spacers[1])
  b2end   = b2start + bclength
  b1len = b1end - b1start

  # filtering and logging
  if filterShortandLongBarcodes(b1len,b2end,bcseq,counts) == 'fail': return None
  logExactOrRegexMatch(spacers,oligo,counts)
  logFuzzyMatching(b1len,bclength,spacers,oligo,counts)

  return [b1start,b1end,b2start,b2end]

def get_qual_scores(qualstring):
    """ Convert FASTQ quality scores to their integer Q score equivalent """
    quality_list = [ord(x) - 33 for x in qualstring]
    return quality_list

def get_err_prob(Q):
    """ Returns the probability of a given base call being incorrect, based on quality score """
    return(10**(-Q/10))

def barcode_quality_check(qualstring, parameters):
    """
    Input: string representing barcode quality and quality check parameters
    Output: 0 if barcode fails quality check
            1 if barcode passes quality check
    """
    quality_list = get_qual_scores(qualstring)
    number_below_min = sum([x < parameters[0] for x in quality_list])
    average_quality = sum(quality_list)/len(quality_list)
    if number_below_min > parameters[1] or average_quality < parameters[2]:
        return 0
    else:
        return 1

def find_most_common_dcrseq2(dcretclist):
    """
    Counts occurrences of dcr-seq pairs in dcretc_list. 
    If one occurs the most, returns this dcr-seq pair
    If more than one occurs the maximum number of times, returns the one with the highest avg seq quality
    If more than one of the maximum occurring dcr-seqs has the same highest avg seq quality,
        returns one of these chosen at random
    """
    dcrseq_count = coll.Counter()
    for detc in dcretclist:
        d = detc['dcr']
        s = detc['seq']
        dcrseq_count['|'.join([d, s])] += 1

    max_value = max(dcrseq_count.values())
    most_common_dcrseqs = [ds for ds, c in dcrseq_count.items() if c == max_value]
    if len(most_common_dcrseqs) == 1:
        return most_common_dcrseqs[0]

    max_avg_quals = []
    for ds in most_common_dcrseqs:
        seq_qual_strings = [detc['seq_qualstring'] for detc in dcretclist if '|'.join([detc['dcr'],detc['seq']]) == ds]
#        seq_qual_strings = [detc['seq_qualstring'] for detc in dcretclist if '|'.join(detc.split('|')[0:2]) == ds]
        seq_avg_quals = [sum(get_qual_scores(x))/len(x) for x in seq_qual_strings]
        max_avg_quals.append(max(seq_avg_quals))
    
    max_indices = [i for i, x in enumerate(max_avg_quals) if x == max(max_avg_quals)]

    if len(max_indices) == 1:
        return most_common_dcrseqs[max_indices[0]]
    
    return most_common_dcrseqs[random.choice(max_indices)]





def find_most_common_dcrseq(dcretclist):
    """
    Counts occurrences of dcr-seq pairs in dcretc_list. 
    If one occurs the most, returns this dcr-seq pair
    If more than one occurs the maximum number of times, returns the one with the highest avg seq quality
    If more than one of the maximum occurring dcr-seqs has the same highest avg seq quality,
        returns one of these chosen at random
    """
    dcrseq_count = coll.Counter()
    for detc in dcretclist:
        d = detc.split('|')[0]
        s = detc.split('|')[1]
        dcrseq_count['|'.join([d, s])] += 1

    max_value = max(dcrseq_count.values())
    most_common_dcrseqs = [ds for ds, c in dcrseq_count.items() if c == max_value]
    if len(most_common_dcrseqs) == 1:
        return most_common_dcrseqs[0]

    max_avg_quals = []
    for ds in most_common_dcrseqs:
        seq_qual_strings = [detc.split('|')[2] for detc in dcretclist if '|'.join(detc.split('|')[0:2]) == ds]
        seq_avg_quals = [sum(get_qual_scores(x))/len(x) for x in seq_qual_strings]
        max_avg_quals.append(max(seq_avg_quals))
    
    max_indices = [i for i, x in enumerate(max_avg_quals) if x == max(max_avg_quals)]
    
    if len(max_indices) == 1:
        return most_common_dcrseqs[max_indices[0]]
    
    return most_common_dcrseqs[random.choice(max_indices)]

def are_seqs_equivalent(seq1, seq2, levdistance_percent_threshold):
    """
    Returns 1 if dcretcs can be considered the same, 0 otherwise
    Defn of equivalent:
        - if sequences are identical OR
        - for identical lengthed sequences, if levenshtein distance as a percentage of length is < threshold
        - for different lengthed sequences, if trimming (from either side) gives a %age lev distance < threshold
    """
    if seq1 == seq2:
        return 1
         
    if len(seq1) == len(seq2):
        lev_percent = 100 * lev.distance(seq1, seq2)/len(seq1)
        if lev_percent <= levdistance_percent_threshold:
            return 1
    
    if len(seq1) != len(seq2):
        lev_percent = 100 * lev.distance(seq1, seq2)/min(len(seq1), len(seq2))
        if lev_percent <= levdistance_percent_threshold:
            return 1

        max_seq_len = max(len(seq1), len(seq2))
        min_seq_len = min(len(seq1), len(seq2))
        both_seqs = [seq1, seq2]
        longerseq = [x for x in both_seqs if len(x) == max_seq_len][0]
        shorterseq = [x for x in both_seqs if len(x) == min_seq_len][0]
        len_diff = len(longerseq) - len(shorterseq)
        
        longerseq_trim1 = longerseq[:len(shorterseq)]
        lev_percent_trim1 = 100 * lev.distance(longerseq_trim1, shorterseq)/len(shorterseq)
        if lev_percent_trim1 <= levdistance_percent_threshold:
            return 1
        
        longerseq_trim2 = longerseq[len_diff:]
        lev_percent_trim2 = 100 * lev.distance(longerseq_trim2, shorterseq)/len(shorterseq)
        if lev_percent_trim2 <= levdistance_percent_threshold:
            return 1
        
    return 0

def h_Dist_Prob():
  # calculate probabilites that a certain Hamming Distance will occur for a comparison between two random UMIs
  percomparison = []
  for HD in range(0,12):
    binomial = comb(12,HD) * 0.75**HD * 0.25**(12-HD)
    percomparison.append(binomial)
  return percomparison

def calc_HD_threshold(n_UMIs, percomparison):
    ncomp = (n_UMIs**2-n_UMIs)/2 # total number of comparisons: half matrix - diagonal
    totalprob = 0.0
    for HD in range(0,12):
      totalprob += ncomp * percomparison[HD] # probability of encountering a certain HD is given by the total number of comparisons times the probability for one comparison
      if totalprob > 0.05: # check if current HD is expected (p>0.05)
        return max(0, HD - 1) # current HD is expected, so return previous HD as threshold, but minimum th=0
    
    print("Error: Could not calculate hamming distance thresholds") # something must have gone wrong if no value is returned at this point
    sys.exit()

def cluster_dcr_barcodes(counter, percomparison, HD_th):
  """
  Takes a counter {barcode:x, barcode:x, ...} of the number of copies of each barcode associated with a dcr
  Returns a counter {barcode_cluster:number_copies, ...} after combining barcodes that are below a threshold
  """
  # make a copy of counter to merge UMIs
  barcodecluster_count = copy.deepcopy(counter)
  if len(counter) == 1:
    # only one barcode is associated with this dcr, no clustering needed
    return counter

  n_UMIs = len(counter) # current number of UMIs during cleaning
  track_barcodes_deleted = coll.Counter()    

  for curr_HD in range(1,12): # go up from HD=1
    n_UMIs = len(barcodecluster_count) # check current number of UMIs

    if n_UMIs == 1:
      return barcodecluster_count

    if HD_th.get(n_UMIs) == None:
      # calculate threshold if not calculated before for this number of clones
      curr_th = calc_HD_threshold(n_UMIs, percomparison)
      HD_th[n_UMIs] = curr_th
    else:
      curr_th = HD_th[n_UMIs]

    if curr_HD > curr_th:
      return barcodecluster_count # if UMIs on this distance are expected to occur, stop merging

    track_barcodes_processed = coll.Counter()

    for bc1, size1 in counter.most_common(): # go from most common to less common UMIs
      if track_barcodes_deleted[bc1] == 0: # check if UMI is not already deleted
        track_barcodes_processed[bc1] = 1 # prevent that UMIs are treated twice during one iteration
        for bc2, size2 in counter.most_common():
          if track_barcodes_deleted[bc2] == 0 and track_barcodes_processed[bc2] == 0:
            if are_barcodes_equivalent(bc1, bc2, curr_HD):
              track_barcodes_deleted[bc2] = 1 # mark smaller UMI as deleted
              track_barcodes_processed[bc2] = 1 # mark smaller UMI as processed
              barcodecluster_count[bc1] += size2 # add size of smaller UMI to size of larger UMI
              del barcodecluster_count[bc2] # remove erroneous UMI from counter
 
  print("Error: Barcode clustering failed")  # something must have gone wrong if no value is returned at this point
  sys.exit()

def NEWcluster_dcr_barcodes(cluster, percomparison, HD_th):
  """
  Takes a counter {barcode:x, barcode:x, ...} of the number of copies of each barcode associated with a dcr
  Returns a counter {barcode_cluster:number_copies, ...} after combining barcodes that are below a threshold
  """
  # make a copy of counter to merge UMIs

  counter = coll.Counter()
  for m in cluster["members"]: 
    counter[m["barcode"]] += 1

  barcodecluster_count = copy.deepcopy(counter)
 # barcodecluster_count = len(cluster["members"])
  if len(counter) == 1:
  #if barcodecluster_count == 1:
    # only one barcode is associated with this dcr, no clustering needed
    return counter
  #return barcodecluster_count

  n_UMIs = len(counter) # current number of UMIs during cleaning
#  n_UMIs = barcodecluster_count # current number of UMIs during cleaning
  track_barcodes_deleted = coll.Counter()    

  for curr_HD in range(1,12): # go up from HD=1
    n_UMIs = len(barcodecluster_count) # check current number of UMIs
    #n_UMIs = barcodecluster_count # check current number of UMIs

    if n_UMIs == 1:
      return barcodecluster_count

    if HD_th.get(n_UMIs) == None:
      # calculate threshold if not calculated before for this number of clones
      curr_th = calc_HD_threshold(n_UMIs, percomparison)
      HD_th[n_UMIs] = curr_th
    else:
      curr_th = HD_th[n_UMIs]

    if curr_HD > curr_th:
      return barcodecluster_count # if UMIs on this distance are expected to occur, stop merging

    track_barcodes_processed = coll.Counter()

    for bc1, size1 in counter.most_common(): # go from most common to less common UMIs
      if track_barcodes_deleted[bc1] == 0: # check if UMI is not already deleted
        track_barcodes_processed[bc1] = 1 # prevent that UMIs are treated twice during one iteration
        for bc2, size2 in counter.most_common():
          if track_barcodes_deleted[bc2] == 0 and track_barcodes_processed[bc2] == 0:
            if are_barcodes_equivalent(bc1, bc2, curr_HD):
              track_barcodes_deleted[bc2] = 1 # mark smaller UMI as deleted
              track_barcodes_processed[bc2] = 1 # mark smaller UMI as processed
              barcodecluster_count[bc1] += size2 # add size of smaller UMI to size of larger UMI
              del barcodecluster_count[bc2] # remove erroneous UMI from counter
 
  print("Error: Barcode clustering failed")  # something must have gone wrong if no value is returned at this point
  sys.exit()


def are_barcodes_equivalent(bc1, bc2, threshold):
    if lev.distance(bc1, bc2) <= threshold:
        return 1
    else:
        return 0

def makecounter():
    return coll.Counter()

def read_in_data_OLD(barcode_quality_parameters, infile):
    ###########################################
    ############# READING DATA IN #############
    ###########################################        
        
    # Check whether file appears to contain suitable verbose Decombinator output for collapsing
    if inputargs['dontcheckinput'] == False:
      if check_dcr_file(infile) != True:
        print("Please check that file contains suitable Decombinator output for collapsing.")
        print("Alternatively, disable the input file sanity check by changing the \'dontcheckinput\' flag, i.e. \'-di True\'")
        sys.exit()
       
    print('Reading data in...')
    t0 = time()
    barcode_dcretc = coll.defaultdict(list)
    #dcr_data = []

    input_dcr_counts = coll.Counter()
    inhandle = opener(infile, 'rt')
    
    for line in inhandle:
        counts['readdata_input_dcrs'] += 1
        fields = line.rstrip('\n').split(', ')

        bc_locs = get_barcode(fields[8],inputargs,counts)        # barcode locations

        if bc_locs:
            # account for N1 barcode being greater or shorter than 6 nt (due to manufacturing errors)
            if (bc_locs[1] - bc_locs[0]) == 6:
                barcode = fields[8][bc_locs[0]:bc_locs[1]] + fields[8][bc_locs[2]:bc_locs[3]]
                barcode_qualstring = fields[9][bc_locs[0]:bc_locs[1]] + fields[9][bc_locs[2]:bc_locs[3]]
            
            elif (bc_locs[1] - bc_locs[0]) < 6:
                n1_diff_len = 6 - (bc_locs[1] - bc_locs[0])
                barcode = fields[8][bc_locs[0]:bc_locs[1]] + "S" * n1_diff_len + fields[8][bc_locs[2]:bc_locs[3]]
                barcode_qualstring = fields[9][bc_locs[0]:bc_locs[1]] + "?" * n1_diff_len + fields[9][bc_locs[2]:bc_locs[3]]
                
                counts['readdata_short_barcode'] += 1

            elif (bc_locs[1] - bc_locs[0]) > 6:
                n1_diff_len = 6 - (bc_locs[1] - bc_locs[0])
                barcode = fields[8][bc_locs[0]:bc_locs[0]+5] + "L" + fields[8][bc_locs[2]:bc_locs[3]]
                barcode_qualstring = fields[9][bc_locs[0]:bc_locs[0]+5] + "?" * n1_diff_len + fields[9][bc_locs[2]:bc_locs[3]]
                
                counts['readdata_long_barcode'] += 1          

                # L and S characters get quality scores of "?", representative of Q30 scores

            if not barcode_quality_check(barcode_qualstring, barcode_quality_parameters):
                # barcode is not sufficient quality, skip to next line of file
                counts['readdata_fail_low_barcode_quality'] += 1
                continue
              
            dcr = ', '.join(line.rstrip('\n').split(', ')[:5])
            input_dcr_counts[str(dcr)] += 1
            
            seq = fields[6]
            
            if len(seq) > inputargs['lenthreshold']: 
                # end V tag to start J tag too long to be sane
                counts['readdata_fail_overlong_intertag_seq'] += 1
                continue
            
            counts['readdata_success'] += 1
            seq_qualstring = fields[7]
            seq_id = fields[5]
            dcretc = '|'.join([dcr, seq, seq_qualstring, seq_id])
            #dcr_data.append(dcretc)
            barcode_dcretc[barcode].append(dcretc)
        else:
            counts['readdata_fail_no_bclocs'] += 1

    inhandle.close()
    
  #  counts['readdata_barcode_dcretc_keys'] = len(barcode_dcretc.keys())
    counts['number_input_unique_dcrs'] = len(input_dcr_counts.keys())
    counts['number_input_total_dcrs'] = sum(input_dcr_counts.values())
      
    t1 = time()
    print('  ', round(t1-t0, 2), 'seconds')
    counts['time_readdata_s'] = t1

    return barcode_dcretc

def read_in_data(barcode_quality_parameters, infile, lev_threshold):
    ###########################################
    ############# READING DATA IN #############
    ###########################################        
        
    # Check whether file appears to contain suitable verbose Decombinator output for collapsing
    if inputargs['dontcheckinput'] == False:
      if check_dcr_file(infile) != True:
        print("Please check that file contains suitable Decombinator output for collapsing.")
        print("Alternatively, disable the input file sanity check by changing the \'dontcheckinput\' flag, i.e. \'-di True\'")
        sys.exit()
       
    print('Reading data in...')
    t0 = time()
    # barcode_dcretc = coll.defaultdict(list)
    # input data will initially be grouped by barcode
#    data = coll.defaultdict(list)
    data = []
    data_barcode_lookup = coll.defaultdict(list)

    input_dcr_counts = coll.Counter()
    inhandle = opener(infile, 'rt')
    
    for line in inhandle:
        counts['readdata_input_dcrs'] += 1
        fields = line.rstrip('\n').split(', ')

        bc_locs = get_barcode(fields[8],inputargs,counts)        # barcode locations

        if bc_locs:
            # account for N1 barcode being greater or shorter than 6 nt (due to manufacturing errors)
            if (bc_locs[1] - bc_locs[0]) == 6:
                barcode = fields[8][bc_locs[0]:bc_locs[1]] + fields[8][bc_locs[2]:bc_locs[3]]
                barcode_qualstring = fields[9][bc_locs[0]:bc_locs[1]] + fields[9][bc_locs[2]:bc_locs[3]]
            
            elif (bc_locs[1] - bc_locs[0]) < 6:
                n1_diff_len = 6 - (bc_locs[1] - bc_locs[0])
                barcode = fields[8][bc_locs[0]:bc_locs[1]] + "S" * n1_diff_len + fields[8][bc_locs[2]:bc_locs[3]]
                barcode_qualstring = fields[9][bc_locs[0]:bc_locs[1]] + "?" * n1_diff_len + fields[9][bc_locs[2]:bc_locs[3]]
                
                counts['readdata_short_barcode'] += 1

            elif (bc_locs[1] - bc_locs[0]) > 6:
                n1_diff_len = 6 - (bc_locs[1] - bc_locs[0])
                barcode = fields[8][bc_locs[0]:bc_locs[0]+5] + "L" + fields[8][bc_locs[2]:bc_locs[3]]
                barcode_qualstring = fields[9][bc_locs[0]:bc_locs[0]+5] + "?" * n1_diff_len + fields[9][bc_locs[2]:bc_locs[3]]
                
                counts['readdata_long_barcode'] += 1          

                # L and S characters get quality scores of "?", representative of Q30 scores

            if not barcode_quality_check(barcode_qualstring, barcode_quality_parameters):
                # barcode is not sufficient quality, skip to next line of file
                counts['readdata_fail_low_barcode_quality'] += 1
                continue
              
            dcr = ', '.join(line.rstrip('\n').split(', ')[:5])
            input_dcr_counts[str(dcr)] += 1
            
            seq = fields[6]
            
            if len(seq) > inputargs['lenthreshold']: 
                # end V tag to start J tag too long to be sane
                counts['readdata_fail_overlong_intertag_seq'] += 1
                continue
            
            counts['readdata_success'] += 1
            seq_qualstring = fields[7]
            seq_id = fields[5]
            dcretc = {'barcode': barcode, 'dcr': dcr, 'seq': seq,
                      'seq_qualstring': seq_qualstring, 'seq_id': seq_id,
                      'collapsed': False}
          
            group_assigned = False
            
            if barcode in data_barcode_lookup:

              for index in data_barcode_lookup[barcode]:
                group = data[index]
                if are_seqs_equivalent(group['protoseq'], dcretc['seq'], lev_threshold):
                  group['members'].append(dcretc)       
                  protodcrseq = find_most_common_dcrseq2(group['members'])
                  group['protodcr'] = protodcrseq.split('|')[0]
                  group['protoseq'] = protodcrseq.split('|')[1]
                  group_assigned = True
                  break

            if not group_assigned:
              data_group = {'barcode': barcode, 'protodcr': dcr, 'protoseq': seq, 'clustered' : False, 'members' : [dcretc]}          
              data.append(data_group)
              data_barcode_lookup[barcode].append(len(data)-1)
              group_assigned = True

        else:
            counts['readdata_fail_no_bclocs'] += 1

    inhandle.close()
    
    counts['number_data_groups'] = len(data)
    counts['number_input_unique_dcrs'] = len(input_dcr_counts.keys())
    counts['number_input_total_dcrs'] = sum(input_dcr_counts.values())
      
    t1 = time()
    print('  ', round(t1-t0, 2), 'seconds')
    counts['time_readdata_s'] = t1

    return data


def cluster_UMIs(data):
    # input data of form: {barcode: [{member1}, {member2},...], barcode: [{member1}, {member2},...]},
    # where each member is a dictionary based on line from input file
    # A member has structure {'barcode': str, 'dcr': str, 'seq': str, 'seq_qualstring': str, 'seq_id': str, 'collapsed': bool}

    total_data = len(data)
    clusters = []

    print("Clustering barcodes...")
    t0 = time()

    for i, group in enumerate(data):
    
      if i%1000 == 0:
        print(i,"/",total_data, "( time elapsed:", round(time() - t0,2),")")

      for j, cluster in enumerate(clusters):  
      
        if not are_barcodes_equivalent(group["barcode"], cluster["barcode"], 2):
          continue

        if not are_seqs_equivalent(group["protoseq"], cluster["protoseq"], 8):
          continue

        cluster["members"] += group["members"]
        group["clustered"] = True
        protodcrseq = find_most_common_dcrseq2(cluster['members'])
        cluster['protodcr'] = protodcrseq.split('|')[0]
        cluster['protoseq'] = protodcrseq.split('|')[1]

        break

      if not group["clustered"]:
        clusters.append(group)
        group["clustered"] = True

    t1 = time()
    print('  ', round(t1-t0, 2), 'seconds')   

    ###############
    # Analysis for testing
    # print("")
    # print("Clusters with over 50 members:")
    # over50 = 0
    # for i,k in enumerate(clusters): 
    #   if len(k["members"]) > 1: 
    #     print (i, len(k["members"]) ) 
    #     over50 += 1
    # if over50 == 0:
    #   print("\tNone found")

    # seq_ids = {} 
    # for i,k in enumerate(clusters): 
    #   seq_ids[i] = set(map(lambda x:x["seq_id"].split("|")[0], k["members"]))


    # print("")
    # print("Clusters with under 10 members:")
    # for k, clus in enumerate(clusters): 
    #   if len(clus["members"]) < 10: 
    #     # print (k, len(clus["members"]) ) 
    #     idx = clusters[k]["members"][0]["seq_id"].split("|")[0] 
    #     #print("idx:",idx) 
    #     clusters_found_in = []
    #     for s in seq_ids:  
    #       if idx in seq_ids[s]:  
    #         clusters_found_in.append(s)

    #     if len(clusters_found_in) > 1:
    #       print (k, len(clus["members"]) ) 
    #       print("idx:",idx) 
    #       for c in clusters_found_in:
    #         found = clusters[c]
    #         bcdist = lev.distance(found["barcode"], clus["barcode"])
    #         seq_in_10 = are_seqs_equivalent(found["protoseq"], clus["protoseq"], 10)
    #         seq_in_5 = are_seqs_equivalent(found["protoseq"], clus["protoseq"], 5)
    #         seq_in_2 = are_seqs_equivalent(found["protoseq"], clus["protoseq"], 2)

    #         print(c, "barcode dist:", bcdist, "seqin10:", seq_in_10, "seqin5:", seq_in_5, "seqin2:", seq_in_2)
    #       print("")
    #############

    # clusters = []
    # percomparison = h_Dist_Prob()

    # print("Clustering barcodes...")
    # t0 = time()

    # # mem_total = 0
    # # for d in data:
    # #   mem_total += len(d["members"])

    # total = len(data)

    # for i, group1 in enumerate(data):

    #   if i%100 == 0:
    #     print("clustered", str(i) + "/" + str(total), "time:", str(time() - t0) )

    #   if group1["clustered"]:
    #     continue
    #   # c_mem_sum = 0
    #   # for c in clusters:
    #   #   c_mem_sum += len(c["members"])
    #   # print("Total members:", mem_total)
    #   # print("No. of members in clusters:", c_mem_sum)

    #   candidates = []

    #   for j, group2 in enumerate(data[i+1:]):

    #     # if i == j:
    #     #   continue

    #     if group2["clustered"]:
    #       continue

    #     # remove non-matching barcodes (less than max threshold) early, for speed
    #     # if not are_barcodes_equivalent(group1["barcode"], group2["barcode"], 7):
    #     #   continue

    #     # seq1 = group1["protodcr"].split(",")[4]
    #     # seq2 = group2["protodcr"].split(",")[4] 

    #     # equiv_seq = are_seqs_equivalent(seq1, seq2, 5)

    #     equiv_seq = are_barcodes_equivalent(group1["protoseq"], group2["protoseq"], 15)

    #     if equiv_seq:
    #       candidates.append(group2)

    #   n_UMIs = len(candidates) + 1     # candidates plus group1

    #   if n_UMIs == 1:

    #     clusters.append(group1)
    #     group1["clustered"] = True
    #     continue

    #   bc_threshold = calc_HD_threshold(n_UMIs, percomparison)
    #   # print("n_UMIS:", n_UMIs)
    #   # print("bc_thresh:", bc_threshold)
    #   # print("")

    #   for k, candidate in enumerate(candidates):

    #     equiv_barcode = are_barcodes_equivalent(group1["barcode"], candidate["barcode"], bc_threshold)

    #     if equiv_barcode:

    #       group1["members"] += candidate["members"]
    #       candidate["clustered"] = True
    #       protodcrseq = find_most_common_dcrseq2(group1['members'])
    #       group1['protodcr'] = protodcrseq.split('|')[0]
    #       group1['protoseq'] = protodcrseq.split('|')[1]


    #   clusters.append(group1)
    #   group1["clustered"] = True

    # t1 = time()
    # print('  ', round(t1-t0, 2), 'seconds')

    # embed()

    # total_data = len(data)

    # for i,group in enumerate(data):

    #   if i%100 == 0:
    #     print(i,"/",total_data)

    #   for cluster in clusters:

    #     equiv_seq = are_seqs_equivalent(group["protoseq"], cluster["protoseq"], lev_threshold)

    #     if not equiv_seq:
    #       continue

    #     # n_UMIs = max(len(group["members"]), len(cluster["members"]))
    #     # percomparison = h_Dist_Prob()

    #     # if n_UMIs == 1:
    #     #   embed()

    #     # bc_threshold = calc_HD_threshold(n_UMIs, percomparison)

        

    #     equiv_barcode = are_barcodes_equivalent(group["barcode"], cluster["barcode"], 5)
     
    #     if equiv_barcode and equiv_seq:

    #       cluster["members"] += group["members"]
    #       # recompute most common dcrseq due to added members
    #       protodcrseq = find_most_common_dcrseq2(cluster['members'])
    #       cluster['protodcr'] = protodcrseq.split('|')[0]
    #       cluster['protoseq'] = protodcrseq.split('|')[1]
    #       group["clustered"] = True
    #       break

    #   if not group["clustered"]:
    #     clusters.append(group)
    #     group["clustered"] = True

    # t1 = time()
    # print('  ', round(t1-t0, 2), 'seconds')

    return clusters

def collapsinate(barcode_quality_parameters,
                 lev_threshold,
                 barcode_distance_threshold,
                 infile, 
                 outpath,
                 file_id):

  
    # read in, structure, and quality check input data
    data = read_in_data(barcode_quality_parameters, infile, lev_threshold)

    # cluster similar UMIs

    clusters  = cluster_UMIs(data)

    # collapse (count) UMIs in each cluster and print to output file

    counts['number_output_unique_dcrs'] = len(clusters)
    counts['number_output_total_dcrs'] = sum(map(lambda x: len(x["members"]),clusters))

    
    # for dcr, copies in dcr_originalcopies.items():
    #     print(', '.join([dcr, str(copies)]), file=outhandle)
    # outhandle.close()

    print("Collapsing clusters...")
    t0 = time()

    collapsed = coll.Counter()
    # true_ids = coll.defaultdict(list)

    for cluster in clusters:
      collapsed[cluster["protodcr"]] += 1
      # seq_ids = set(map(lambda x: x['seq_id'],clusters[111]["members"]))
      # true_ids[cluster["protodcr"]] += list(seq_ids)
      # true_ids[cluster["protodcr"]] = list(set(true_ids[cluster["protodcr"]])) 
    
    t1 = time()
    print('  ', round(t1-t0, 2), 'seconds')  

    outfile = outpath + file_id + suffix
    outhandle = open(outfile, 'w')
    print('Writing to output file', outfile, '...')

    for dcr, dcr_count in collapsed.items():
      print(', '.join([dcr, str(dcr_count)]), file=outhandle)
    outhandle.close()

    if inputargs['dontgzip'] == False:  # Gzip output file
        print("Compressing output to", outfile + ".gz ...")
      
        with open(outfile) as inf, gzip.open(outfile + '.gz', 'wt') as outf:
            outf.writelines(inf)
        outf.close()
        os.unlink(outfile)

        outfilenam = outfile + ".gz"
    else:
        outfilenam = outfile + suffix    

    counts['outfile'] = outfilenam

    # only need to run this bit if interested in the number of times each barcode is repeated in the data
    if inputargs['barcodeduplication'] == True:
        outfile = outpath+file_id+'_barcode_duplication.txt'
        outhandle = open(outfile, 'w')
        for bc, copies in barcode_duplication.items():
            print(','.join([bc, str(copies)]), file=outhandle)
        outhandle.close()
        
    counts['outfilenam'] = outfilenam

    return 1



def collapsinate_SEMIOLD(barcode_quality_parameters,
                 lev_threshold,
                 barcode_distance_threshold,
                 infile, 
                 outpath,
                 file_id):

    # read in, structure, and quality check input data
    dcr_clusters = read_in_data(barcode_quality_parameters, infile)

    ##########################################
    ############# TCR COLLAPSING #############
    ##########################################
    
    print('Collapsing by barcode...')
    t0 = time()

    clusters = []

    for barcode, dcrs in dcr_clusters.items():
      
      # from IPython import embed
      # embed()
      number_collapsed = 0
     # protoseqcount = 0
      while number_collapsed < len(dcrs):
        uncollapsed_dcrs = [i for i in filter(lambda x: x["collapsed"] != True, dcrs)] 
        protodcrseq = find_most_common_dcrseq2(uncollapsed_dcrs)
        protodcr = protodcrseq.split('|')[0]
        protoseq = protodcrseq.split('|')[1] 
        #protoseqcount += 1

        cluster = {"protodcr": protodcr, "protoseq": protoseq,"barcode": barcode, "members": [], "combined": False}

        for di, dcretc in enumerate(uncollapsed_dcrs):
          if are_seqs_equivalent(protoseq, dcretc['seq'], lev_threshold):
            cluster["members"].append(dcretc)
            dcretc["collapsed"] = True
            number_collapsed += 1
          ## NEW FUNCTIONALITY : SAME V AND J, SIMILAR INSERT
          else:
            # dcretc gets left in remaining_dcretc_list, and is considered next time around
            # print(protoseq)
            # print(dcretc['seq'])
            # from IPython import embed
            # embed()
            pass

        # if protodcr not in clusters:
        #   clusters[protodcrseq] = [cluster]
        # else:
        #   clusters[protodcrseqs].append(cluster)
        clusters.append(cluster)

    # number of error-corrected TCRs remaining
   # counts['tcrcoll_dcr_barcodecounter_keys'] = len(clusters.keys()) 
   
    t1 = time()
    print('  ', round(t1-t0, 2), 'seconds')
    counts['time_tcrcollapsing_s'] = t1

    # from IPython import embed
    # embed()

    ##############################################
    ############# COMBINING CLUSTERS #############
    ##############################################   

    print("Combining collapsed clusters...")
    t0 = time()

    combined_clusters = []

 #   cluster_threshold = 50
 #   clusters.sort(key = lambda x: len(x["members"]), reverse = True)

    for i, c1 in enumerate(clusters):

      if i%500 == 0:
        print(i)

      # if len(c1["members"]) > cluster_threshold:
      #   combined_clusters.append(c1)
      #   c1["combined"] = True
      #   continue

      # if i == 0:
      #   combined_clusters.append(c1)
      #   c1["combined"] = True

      codes = list(map(lambda x: x["barcode"], combined_clusters))
      close_umi_cluster_indexes = [i for i, bcode in enumerate(codes) if lev.distance(bcode,c1["barcode"]) <= 2]


      for j, cl_index in enumerate(close_umi_cluster_indexes):

        c2 = combined_clusters[cl_index]
        # combine clusters if barcodes are within threshold of levenshtein distance of 1
        # AND protoseqs of clusters are within similarity threshold
        # AND same V and J are used
        # AND insert seqs are within threshold of levenshtein distance of 1
       # eq_barcode = are_barcodes_equivalent(c1["barcode"],c2["barcode"],1)
        eq_seq = are_seqs_equivalent(c1["protoseq"],c2["protoseq"], lev_threshold)
        vj1 = c1["protodcr"].replace(" ","").split(",")[0:2]
        vj2 = c2["protodcr"].replace(" ","").split(",")[0:2]
    #    ins1 = c1["protodcr"].replace(" ","").split(",")[4]
    #    ins2 = c2["protodcr"].replace(" ","").split(",")[4]
        if vj1 == vj2:
          eq_vj = 1
        else:
          eq_vj = 0

        # use barcode equivalency... just a simple levenshtein comparison (rename?)
    #    eq_insert = are_barcodes_equivalent(ins1, ins2, 1)

        if eq_seq and eq_vj:

          # add c2 cluster to c1 cluster in combined clusters
          c2["members"] += c1["members"]
          c1["combined"] = True

          # set new protodcr and protoseq in case they have changed
          protodcrseq = find_most_common_dcrseq2(c2["members"])
          c2["protodcr"] = protodcrseq.split("|")[0]
          c2["protoseq"] = protodcrseq.split("|")[1]
          break

      if not c1["combined"]:
        combined_clusters.append(c1)
        c1["combined"] = True
  
     
    t1 = time()
    print('  ', round(t1-t0, 2), 'seconds')
    
    counts['time_combiningclusters_s'] = t1

    # from IPython import embed
    # embed()
    
    ##############################################
    ############# BARCODE CLUSTERING #############
    ##############################################   

    print('Clustering barcodes...')
    t0 = time()
    
    dcr_originalcopies = coll.Counter()         # Stores the final collapsed DCRs with estimated frequencies (i.e. # clustered barcodes)
    barcode_duplication = coll.Counter()        # Stores the frequency with which each barcode appeared in input file (number lines)
    
    # calculate for a given number of UMIs the maximum Hamming distance that is not expected to occur (p<0.05)
    h_probs = h_Dist_Prob()
    HD_th = {}

    # for dcr, barcode_count in dcr_barcodecounter.items():
    for clus in combined_clusters:

        dcr = clus["protodcr"]
        #clustered_bc_dict = OLDcluster_dcr_barcodes(barcode_count, barcode_distance_threshold)
      #  clustered_bc_dict = cluster_dcr_barcodes(barcode_count, h_probs, HD_th)
        clustered_bc_dict = NEWcluster_dcr_barcodes(clus, h_probs, HD_th)

        dcr_originalcopies[dcr] = len(clustered_bc_dict)
        for bc, c in clustered_bc_dict.items():
            barcode_duplication[bc] += c
    
    counts['number_output_unique_dcrs'] = len(dcr_originalcopies)
    counts['number_output_total_dcrs'] = sum(dcr_originalcopies.values())
    
    t1 = time()
    counts['time_bcclustering_s'] = t1
    print('  ', round(t1-t0, 2), 'seconds')

    # from IPython import embed
    # embed()

    outfile = outpath + file_id + suffix
    outhandle = open(outfile, 'w')
    print('Writing to output file', outfile, '...')
    
    for dcr, copies in dcr_originalcopies.items():
        print(', '.join([dcr, str(copies)]), file=outhandle)
    outhandle.close()
    
    if inputargs['dontgzip'] == False:  # Gzip output file
        print("Compressing output to", outfile + ".gz ...")
      
        with open(outfile) as inf, gzip.open(outfile + '.gz', 'wt') as outf:
            outf.writelines(inf)
        outf.close()
        os.unlink(outfile)

        outfilenam = outfile + ".gz"
    else:
        outfilenam = outfile + suffix    

    counts['outfile'] = outfilenam

    # only need to run this bit if interested in the number of times each barcode is repeated in the data
    if inputargs['barcodeduplication'] == True:
        outfile = outpath+file_id+'_barcode_duplication.txt'
        outhandle = open(outfile, 'w')
        for bc, copies in barcode_duplication.items():
            print(','.join([bc, str(copies)]), file=outhandle)
        outhandle.close()
        
    counts['outfilenam'] = outfilenam



def collapsinate_OLD(barcode_quality_parameters,
                 lev_threshold,
                 barcode_distance_threshold,
                 infile, 
                 outpath,
                 file_id):
  
    # read in, structure, and quality check input data
    barcode_dcretc = read_in_data_OLD(barcode_quality_parameters, infile)

    ##########################################
    ############# TCR COLLAPSING #############
    ##########################################        
    
    print('Collapsing by barcode...')
    t0 = time()
    
    dcr_barcodecounter = coll.defaultdict(makecounter)

    for barcode, dcretc_list in barcode_dcretc.items():

        remaining_dcretc_list = dcretc_list[:]

        while remaining_dcretc_list:
            protodcrseq = find_most_common_dcrseq(remaining_dcretc_list)
            protodcr = protodcrseq.split('|')[0]
            protoseq = protodcrseq.split('|')[1]
  
            for di, dcretc in enumerate(remaining_dcretc_list):
                thisseq = dcretc.split('|')[1]
                if are_seqs_equivalent(protoseq, thisseq, lev_threshold):
                    dcr_barcodecounter[protodcr][barcode] += 1
                    counts['tcrcoll_kept_bc'] += 1
                    remaining_dcretc_list.pop(di)
                else:
                    # dcretc gets left in remaining_dcretc_list, and is considered next time around
                    pass

    # number of error-corrected TCRs remaining
    counts['tcrcoll_dcr_barcodecounter_keys'] = len(dcr_barcodecounter.keys()) 
     
    t1 = time()
    print('  ', round(t1-t0, 2), 'seconds')
    counts['time_tcrcollapsing_s'] = t1
    
    ##############################################
    ############# BARCODE CLUSTERING #############
    ##############################################          

    print('Clustering barcodes...')
    t0 = time()
    
    dcr_originalcopies = coll.Counter()         # Stores the final collapsed DCRs with estimated frequencies (i.e. # clustered barcodes)
    barcode_duplication = coll.Counter()        # Stores the frequency with which each barcode appeared in input file (number lines)
    
    # calculate for a given number of UMIs the maximum Hamming distance that is not expected to occur (p<0.05)
    h_probs = h_Dist_Prob()
    HD_th = {}

    for dcr, barcode_count in dcr_barcodecounter.items():

        #clustered_bc_dict = OLDcluster_dcr_barcodes(barcode_count, barcode_distance_threshold)
        clustered_bc_dict = cluster_dcr_barcodes(barcode_count, h_probs, HD_th)
        dcr_originalcopies[dcr] = len(clustered_bc_dict)
        for bc, c in clustered_bc_dict.items():
            barcode_duplication[bc] += c
    
    counts['number_output_unique_dcrs'] = len(dcr_originalcopies)
    counts['number_output_total_dcrs'] = sum(dcr_originalcopies.values())
    
    t1 = time()
    counts['time_bcclustering_s'] = t1
    print('  ', round(t1-t0, 2), 'seconds')

    outfile = outpath + file_id + suffix
    outhandle = open(outfile, 'w')
    print('Writing to output file', outfile, '...')
    
    for dcr, copies in dcr_originalcopies.items():
        print(', '.join([dcr, str(copies)]), file=outhandle)
    outhandle.close()
    
    if inputargs['dontgzip'] == False:  # Gzip output file
        print("Compressing output to", outfile + ".gz ...")
      
        with open(outfile) as inf, gzip.open(outfile + '.gz', 'wt') as outf:
            outf.writelines(inf)
        outf.close()
        os.unlink(outfile)

        outfilenam = outfile + ".gz"
    else:
        outfilenam = outfile + suffix    

    counts['outfile'] = outfilenam

    # only need to run this bit if interested in the number of times each barcode is repeated in the data
    if inputargs['barcodeduplication'] == True:
        outfile = outpath+file_id+'_barcode_duplication.txt'
        outhandle = open(outfile, 'w')
        for bc, copies in barcode_duplication.items():
            print(','.join([bc, str(copies)]), file=outhandle)
        outhandle.close()
        
    counts['outfilenam'] = outfilenam
                
########################################################################################################################
if __name__ == '__main__':

    ########################################
    # Get parameters 
    inputargs = vars(args())
    
    if inputargs['infile'].endswith('.gz'):
      opener = gzip.open
    else:
      opener = open    
    
    suffix = "." + inputargs['extension']
    
    counts = coll.Counter()
    counts['start_time'] = time()   
        
    ## this is [min_barcode_nt_quality, max_bc_nts_with_min_quality, min_avg_bc_quality]
    barcode_quality_parameters = [inputargs['minbcQ'], inputargs['bcQbelowmin'], inputargs['avgQthreshold']]

    ## this is the percentage lev distance that is allowed to determine whether two sequences are equivalent
    lev_threshold = inputargs['percentlevdist']
    
    ## this is the number of barcode edits that are allowed to call two barcodes equivalent
    barcode_distance_threshold = inputargs['bcthreshold']

    infile = inputargs['infile']
    outpath = ''
    
    file_id = infile.split('/')[-1].split('.')[0]
    ########################################

    collapsinate(barcode_quality_parameters,
                 lev_threshold,
                 barcode_distance_threshold,
                 infile, 
                 outpath,
                 file_id)
    
    counts['end_time'] = time()    
    counts['time_taken_total_s'] = counts['end_time'] - counts['start_time']
    
    #######################################
    
    # Write data to summary file
    if inputargs['suppresssummary'] == False:
      
      # Check for directory and make summary file
      if not os.path.exists('Logs'):
        os.makedirs('Logs')
      date = strftime("%Y_%m_%d")
      
      # Check for existing date-stamped file
      summaryname = "Logs/" + date + "_" + file_id + "_Collapsing_Summary.csv"
      if not os.path.exists(summaryname): 
        summaryfile = open(summaryname, "w")
      else:
        # If one exists, start an incremental day stamp
        for i in range(2,10000):
          summaryname = "Logs/" + date + "_" + file_id + "_Collapsing_Summary" + str(i) + ".csv"
          if not os.path.exists(summaryname): 
            summaryfile = open(summaryname, "w")
            break
          
      # Generate string to write to summary file

      summstr = "Property,Value\nVersion," + str(__version__) + "\nDirectory," + os.getcwd() + "\nInputFile," + inputargs['infile'] \
        + "\nOutputFile," + counts['outfilenam'] + "\nDateFinished," + date + "\nTimeFinished," + strftime("%H:%M:%S") \
        + "\nTimeTaken(Seconds)," + str(round(counts['time_taken_total_s'],2)) + "\n\n"

      for s in ['extension', 'dontgzip', 'allowNs', 'dontcheckinput', 'barcodeduplication', 'minbcQ', 'bcQbelowmin', 'bcthreshold', \
        'lenthreshold', 'percentlevdist', 'avgQthreshold', 'positionalbarcodes']:
        summstr = summstr + s + "," + str(inputargs[s]) + "\n"

      counts['pc_input_dcrs'] = counts['number_input_total_dcrs'] / counts['readdata_input_dcrs']
      counts['pc_uniq_dcr_kept'] = ( counts['number_output_unique_dcrs'] / counts['number_input_unique_dcrs'] )
      counts['pc_total_dcr_kept'] = ( counts['number_output_total_dcrs'] / counts['number_input_total_dcrs'] )
      
      counts['avg_input_tcr_size'] = counts['number_input_total_dcrs'] / counts['number_input_unique_dcrs']
      counts['avg_output_tcr_size'] = counts['number_output_total_dcrs'] / counts['number_output_unique_dcrs']
      counts['avg_RNA_duplication'] = 1 / counts['pc_total_dcr_kept']
      
      # success properties
      summstr = summstr + "\nInputUncollapsedDCRLines," + str(counts['readdata_input_dcrs']) \
        + "\nUniqueDCRsPassingFilters," + str(counts['number_input_unique_dcrs']) \
        + "\nTotalDCRsPassingFilters," + str(counts['number_input_total_dcrs']) \
        + "\nPercentDCRPassingFilters(withbarcode)," + str( round(counts['pc_input_dcrs'], 3 ) ) \
        + "\nUniqueDCRsPostCollapsing," + str(counts['number_output_unique_dcrs']) \
        + "\nTotalDCRsPostCollapsing," + str(counts['number_output_total_dcrs']) \
        + "\nPercentUniqueDCRsKept," + str( round(counts['pc_uniq_dcr_kept'], 3 ) ) \
        + "\nPercentTotalDCRsKept," + str( round(counts['pc_total_dcr_kept'], 3 ) ) \
        + "\nAverageInputTCRAbundance," + str( round(counts['avg_input_tcr_size'], 3 ) ) \
        + "\nAverageOutputTCRAbundance," + str( round(counts['avg_output_tcr_size'], 3 ) ) \
        + "\nAverageRNAduplication," + str( round(counts['avg_RNA_duplication'], 3 ) ) \
        + "\n\nBarcodeFail_ContainedNs," + str(counts['getbarcode_fail_N']) \
        + "\nBarcodeFail_SpacersNotFound," + str(counts['readdata_fail_no_bclocs']) \
        + "\nBarcodeFail_LowQuality," + str(counts['readdata_fail_low_barcode_quality'])

      print(summstr,file=summaryfile) 
      summaryfile.close()
