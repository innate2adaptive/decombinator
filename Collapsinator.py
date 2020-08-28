"""#######################################################################################################################################################
  Collapsinator
  Katharine Best and James M. Heather, August 2016, UCL
  https://innate2adaptive.github.io/Decombinator/

##################
### BACKGROUND ###
##################

  Takes the output files of Decombinator (run using the barcoding option) and performs collapsing and error correction
  This version is a modified version of KB's script collapsinator_20141126.py
  (That was itself an improved version of the CollapseTCRs.py script used in the Heather et al HIV TCR paper (DOI: 10.3389/fimmu.2015.00644))
  Version 4.0.2 includes improved clustering routines measuring the similarity in both barcode and TCR sequence of TCR repertoire data

##################
###### INPUT #####
##################

  Takes as input .n12 files produced by Decombinator (v3 or higher), assuming it has been run on suitably barcoded and demultiplexed data.

  Other optional flags:
  
    -s/--supresssummary: Supress the production of a summary file containing details of the run into a 'Logs' directory. 
  
    -dz/--dontgzip: Suppress the automatic compression of output demultiplexed FASTQ files with gzip. 
  
    -dc/--dontcount: Suppress the whether or not to show the running line count, every 100,000 reads. Helps in monitoring progress of large batches.

  The other optional flags are somewhat complex, and caution is advised in their alteration.

  To see all options, run: python Collapsinator.py -h

  Input files need to be in the appropriate format, consisting of:
    V index, J index, V deletions, J deletions, insert, ID, inter-tag TCR sequence, inter-tag quality, barcode sequence, barcode quality

##################
##### OUTPUT #####  
##################
   
  A Decombinator index file, giving each error-corrected DCR index, and the frequency with which it appears
  in the final processed data, and an average UMI count, which can be used to estimate the robustness of the
  data for that particular sequence

#######################################################################################################################################################"""

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

__version__ = '4.0.2'

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
      '-dc', '--dontcount', action='store_true', help='Block the line count from being shown', required=False, default=False)
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
  parser.add_argument(
      '-wc', '--writeclusters', action='store_true', help='Write cluster data to separate cluster files',\
        required=False, default=False)
  parser.add_argument(
      '-uh', '--UMIhistogram', action='store_true', help='Creates histogram of average UMI cluster sizes',\
        required=False, default=False)
  
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
  oligos['m13'] = {'spcr1': 'GTCGTGACTGGGAAAACCCTGG','spcr2':'GTCGTGAT'}
  oligos['i8'] = {'spcr1':'GTCGTGAT','spcr2':'GTCGTGAT'}

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
    spcr1 = oligo['spcr1']
    spacer += spacerSearch(spcr1, seq[0:len(spcr1)+allowance])
    return spacer

def findSecondSpacer(oligo,seq):
    spacer = []
    spcr1 = oligo['spcr1']
    spcr2 = oligo['spcr2']
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

def get_barcode_positions(bcseq,inputargs,counts):
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


def set_barcode(fields, bc_locs):
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
    return barcode, barcode_qualstring

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

def are_seqs_equivalent(seq1, seq2, lev_percent_threshold):
  # Returns True if seqs can be considered the same, False otherwise
  # Definition of equivalent:
  #   levenshtein distance as a percentage of the shorter of the two seqs is <= threshold
  threshold = len(min(seq1, seq2, key=len)) * lev_percent_threshold
  return lev.distance(seq1, seq2) <= threshold

def are_barcodes_equivalent(bc1, bc2, threshold):
    return lev.distance(bc1, bc2) <= threshold

def read_in_data(barcode_quality_parameters, infile, lev_threshold, dont_count):
    ###########################################
    ############# READING DATA IN #############
    ###########################################        
        
    # Check whether file appears to contain suitable verbose Decombinator output for collapsing
    if inputargs['dontcheckinput'] == False:
      if check_dcr_file(infile) != True:
        print("Please check that file contains suitable Decombinator output for collapsing.")
        print("Alternatively, disable the input file sanity check by changing the \'dontcheckinput\' flag, i.e. \'-di True\'")
        sys.exit()
       
    print("Reading data in...")
    t0 = time()
    barcode_dcretc = coll.defaultdict(list)
    barcode_lookup = coll.defaultdict(list)

    input_dcr_counts = coll.Counter()
    inhandle = opener(infile, 'rt')
    
    for lcount, line in enumerate(inhandle):
        if lcount % 50000 == 0 and lcount != 0 and not dont_count:
          print("   Read in", lcount, "lines... ", round(time()-t0,2), "seconds")
        counts['readdata_input_dcrs'] += 1
        fields = line.rstrip('\n').split(', ')

        bc_locs = get_barcode_positions(fields[8], inputargs, counts)        # barcode locations

        if not bc_locs:
          counts['readdata_fail_no_bclocs'] += 1
          continue

        barcode, barcode_qualstring = set_barcode(fields, bc_locs)

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
          
        group_assigned = False

        # Assign reads to groups based on their barcode data. Reads with identical barcodes are grouped together
        # so long as they have equivalent TCR sequences. Reads with identical barcodes but non-equivalent TCR
        # sequences are grouped separately.
        # Data is grouped in dictionary format: {'barcode|index|protoseq' : [dcretc1, dcretc2, ...], ... }
        # where index counts upwards from zero to help disinguish identical barcodes in different groups,
        # protoseq is the most common sequence present in the group, and dcretc are the input reads
        if barcode in barcode_lookup:

          for index in barcode_lookup[barcode]:

            if are_seqs_equivalent(index[1], seq, lev_threshold):

              barcode_dcretc[barcode + "|" + str(index[0]) + "|" + index[1]].append(dcretc)       
              protodcretc_list = barcode_dcretc[barcode + "|" + str(index[0]) + "|" + index[1]]
              seq_counter = coll.Counter(map(lambda x: x.split("|")[1],protodcretc_list))
              protoseq = seq_counter.most_common(1)[0][0] # find most common sequence in group
            
              if not index[1] == protoseq:
                # if there is a new protoseq, replace record with old protoseq
                # with identical record with updated  protoseq
                barcode_dcretc[barcode + "|" + str(index[0]) + "|" + protoseq] = barcode_dcretc[barcode + "|" + str(index[0]) + "|" + index[1]]
                del barcode_dcretc[barcode + "|" + str(index[0]) + "|" + index[1]]

                barcode_lookup[barcode][index[0]] = [index[0], protoseq]

              group_assigned = True
              # if assigned to a group, stop and move onto next read
              break

          if not group_assigned:
            # if no appropriate group found, create new group with correctly incremented index
            barcode_lookup[barcode].append([index[0] + 1,seq])
            barcode_dcretc["|".join([barcode,str(index[0]+1),seq])].append(dcretc)
            group_assigned = True

        else:
          # if no identical barcode found, create new barcode group with index zero
          barcode_lookup[barcode].append([0,seq])
          barcode_dcretc["|".join([barcode,"0",seq])].append(dcretc)
          group_assigned = True

    inhandle.close()
    
    counts['readdata_barcode_dcretc_keys'] = len(barcode_dcretc.keys())
    counts['number_input_unique_dcrs'] = len(input_dcr_counts.keys())
    counts['number_input_total_dcrs'] = sum(input_dcr_counts.values())
      
    t1 = time()
    print("   Read in total of", lcount+1, "lines")
    print("  ", counts['readdata_success'], "reads sorted into", len(barcode_dcretc), "initial groups" )
    print('  ', round(t1-t0, 2), 'seconds')
    counts['time_readdata_s'] = t1

    return barcode_dcretc


def cluster_UMIs(barcode_dcretc, inputargs, barcode_threshold, seq_threshold, dont_count):
    # input data of form: {'barcode1|index|protoseq': [dcretc1, dcretc2,...], 'barcode2|index|protoseq|: [dcretc1, dcretc2,...], ...}
    # (see read_in_data function for details)
    # This function merges groups that have both equivalent barcodes and equivalent protoseqs
    # output data of form: {'barcode1|index|protoseq': [dcretc1, dcretc2,...], 'barcode2|index|protoseq|: [dcretc1, dcretc2,...], ...}
    # Output format is same as input format, but after clustering. The protoseq is recalculated when new
    # dcretcs are added to a cluster

    percent_seq_threshold = seq_threshold/100.0

    print("Clustering barcodes groups...")
    t0 = time()

    clusters = coll.defaultdict(list)

    barcode_dcretc_total = len(barcode_dcretc)
    count = 0

    for i, b1 in enumerate(barcode_dcretc):

      clustered = False

      if count % 5000 == 0 and not dont_count:
        print("   Clustered", count, "/", barcode_dcretc_total, "...", round(time()-t0,2),"seconds")

      barcode1, index1, protoseq1 = b1.split("|")
      
      for j,b2 in enumerate(clusters): 
        barcode2, index2, protoseq2 = b2.split("|")

        if are_barcodes_equivalent(barcode1, barcode2, barcode_threshold):

          if are_seqs_equivalent(protoseq1, protoseq2, percent_seq_threshold):

            clusters[b2] += barcode_dcretc[b1]

            protodcretc_list = clusters[b2]
            seq_counter = coll.Counter(map(lambda x: x.split("|")[1],protodcretc_list))
            protoseq = seq_counter.most_common(1)[0][0] # recalculate protoseq (most common sequence), after merging

            if not protoseq2 == protoseq:
              # if there is a new protoseq, replace record with old protoseq
              # with identical record with updated protoseq 
              clusters["|".join([barcode2,index2, protoseq])] = clusters[b2]

              del clusters[b2]

            clustered = True
            break

      if not clustered:
        # if not appropriate cluster found to add group to, then create a new cluster identical to the group
        clusters[b1] = barcode_dcretc[b1]

      count += 1

    t1 = time()
    print("  ", len(barcode_dcretc), "groups merged into", len(clusters), "clusters")
    print("  ", round(t1-t0, 2), "seconds")
    
    # dump clusters to separate files if desired
    if inputargs['writeclusters']:
      write_clusters(clusters)

    return clusters

def write_clusters(clusters):
    # create directory to store cluster data without overwriting exiting directories
    dirname = "clusters"
    count = 1
    while os.path.isdir(dirname):
      dirname = "clusters" + str(count)
      count += 1
    os.mkdir(dirname)

    print("   Writing clusters to directory: ", os.path.abspath(dirname), "...")
    # write data of each cluster to a separate file and store in clusters directory
    for k in clusters: 
      with open(dirname+os.sep+"|".join(k.split("|")[:2])+".txt", 'w') as ofile: 
        for j in clusters[k]: 
          print(j, file = ofile)
    return 1


def collapsinate(barcode_quality_parameters, lev_threshold, barcode_distance_threshold,
                 infile, outpath, file_id, dont_count):
 
    # read in, structure, and quality check input data
    barcode_dcretc = read_in_data(barcode_quality_parameters, infile, lev_threshold, dont_count)

    # cluster similar UMIs
    clusters = cluster_UMIs(barcode_dcretc, inputargs, barcode_distance_threshold, lev_threshold, dont_count)

    # collapse (count) UMIs in each cluster and print to output file
    print("Collapsing clusters...")
    t0 = time()

    collapsed = coll.Counter()
    cluster_sizes = coll.defaultdict(list)

    for c in clusters:
      protodcr = coll.Counter(map(lambda x: x.split("|")[0],clusters[c])).most_common(1)[0][0] # find most common dcr in each cluster
      collapsed[protodcr] += 1
      cluster_sizes[protodcr].append(len(clusters[c]))

    counts['number_output_unique_dcrs'] = len(collapsed)
    counts['number_output_total_dcrs'] = sum(collapsed.values())      

    t1 = time()
    print('  ', round(t1-t0, 2), 'seconds')  

    outfile = outpath + file_id + suffix
    outhandle = open(outfile, 'w')
    print("Writing to output file", outfile, "...")

    average_cluster_size_counter = coll.Counter()
    for dcr, dcr_count in collapsed.items():
      av_clus_size = round(sum(cluster_sizes[dcr])/dcr_count)
      average_cluster_size_counter[av_clus_size] += 1
      print(', '.join([dcr, str(dcr_count), str(av_clus_size)]), file=outhandle)
    outhandle.close()

    if inputargs['dontgzip'] == False:  # Gzip output file
        print("Compressing output to", outfile + ".gz ...")
      
        with open(outfile) as inf, gzip.open(outfile + '.gz', 'wt') as outf:
            outf.writelines(inf)
        outf.close()
        os.unlink(outfile)

        outfilenam = outfile + ".gz"
    else:
        outfilenam = outfile   

    # only need to run this bit if interested in the number of times each barcode is repeated in the data
    if inputargs['barcodeduplication'] == True:
        outfile = outpath+file_id+'_barcode_duplication.txt'
        outhandle = open(outfile, 'w')
        for bc, copies in clusters.items():
            barcode_index = "|".join(bc.split("|")[:2])
            print(','.join([barcode_index, str(len(copies))]), file=outhandle)
        outhandle.close()
        print("barcode duplication data saved to", outfile)
        
    counts['outfilenam'] = outfilenam

    return collapsed, average_cluster_size_counter


if __name__ == '__main__':

    ########################################
    # Get parameters 
    inputargs = vars(args())

    print("Running Collapsinator version", __version__)
    
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

    ## this is a boolean for printing progress of the run to the terminal (False for printing, True for not printing; default = False)
    dont_count = inputargs['dontcount']

    ########################################

    collapsed, average_cluster_size_counter = collapsinate(barcode_quality_parameters,
                                        lev_threshold, barcode_distance_threshold,
                                        infile, outpath, file_id, dont_count)
    
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

      # create output data for UMI histogram and save to file (optional)
      if inputargs['UMIhistogram']:

        hfileprefix = "_".join( summaryname.split("_")[:-2] + ['UMIhistogram'] )
        # iterate to create unique file name
        if os.path.exists(hfileprefix + '.csv'):
          i = 1
          while os.path.exists(hfileprefix + str(i) + '.csv'):
            i += 1
          hfileprefix += str(i)
        
        hfilename = hfileprefix + '.csv'

        with open(hfilename,'w') as hfile:

          for av, count in sorted(average_cluster_size_counter.items()):
            print(str(av) + "," + str(count), file=hfile)

        # print instructions for creating histogram plot using script in Supplementary-Scripts
        print("\nAverage UMI cluster size histogram data saved to", hfilename)
        print("To plot histogram, please use UMIhistogram.py script located in the Decombinator-Tools repository.")
        print("Decombinator-Tools can be found at https://github.com/innate2adaptive/Decombinator-Tools.")
        print("With the Decombinator-Tools repository downloaded, run:")
        codestr = "python path/to/Decombinator-Tools/UMIHistogram.py -in "+ hfilename
        print("#"*(len(codestr) + 4))
        print(" ", codestr, " ")
        print("#"*(len(codestr) + 4))
