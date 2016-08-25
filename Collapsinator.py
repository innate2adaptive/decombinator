# Collapsinator
# Katharine Best and James M. Heather, August 2016, UCL
# https://innate2adaptive.github.io/Decombinator/

##################
### BACKGROUND ###
##################

# Takes the output files of Decombinator (run using the barcoding option) and performs collapsing and error correction
# This version is a modified version of KB's script collapsinator_20141126.py
  # That was itself an improved version of the CollapseTCRs.py script used in the Heather et al HIV TCR paper (DOI: 10.3389/fimmu.2015.00644)

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
import os, sys

__version__ = '2.1'

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
      print 'Cannot find file, please double-check path.'
      return False
    if os.path.getsize(infile) == 0:
      print 'Input file appears to be empty; please double-check path.'
      return False
    
    # Check first few lines    
    with opener(infile) as poss_dcr: 
      for i in range(5):
              
        # Check it's a comma-delimited file
        if "," not in next(poss_dcr):
          print 'Input Decombinator file sanity check fail: seemingly not comma-delimited text file.'
          return False
        else:
          fields = next(poss_dcr).rstrip().split(", ")
        
        # Check lines contain correct number of fields
        if len(fields) != 10:
          print 'Input Decombinator file sanity check fail: file does not contain the correct number of comma-delimited fields (ten).'
          return False
        # Check DCR classifiers are feasible (i.e. 4x integers followed by a DNA string)
        elif not all([num_check(field) for field in fields[0:4]]):
          print 'Input Decombinator file sanity check fail: integer components of Decombinator classifier not feasible (i.e. not integers >= zero).'
          return False        
        elif not is_dna(fields[4]):
          print 'Input Decombinator file sanity check fail: Decombinator insert sequence field contains non-DNA sequence.'
          return False
        # Check inter-tag and barcode sequences are legitimate DNA strings 
        elif is_dna(fields[6]) != True or is_dna(fields[8]) != True:
          print 'Input Decombinator file sanity check fail: inter-tag and/or barcode sequence contains non-DNA sequences.'
          return False
        # Check inter-tag and barcode quality are feasible FASTQ quality scores
        elif all(set([num_check(x) for x in get_qual_scores(fields[7])])) != True or all(set([num_check(x) for x in get_qual_scores(fields[9])])) != True:
          print 'Input Decombinator file sanity check fail: inter-tag and/or barcode quality strings do not appear to contain valid scores.'
          return False        
        # Check that inter-tag and barcode sequence/quality pairs are the correct length
        elif len(fields[6]) != len(fields[7]) or len(fields[8]) != len(fields[9]):
          print 'Input Decombinator file sanity check fail: inter-tag and/or barcode sequence and quality string pairs are not of the same length.'
          return False    
    
      return True       # If first few lines all pass, assume the file is fine and carry on with the analysis.


def get_barcode(bcseq):
    """
    getbarcode(bcseq):
    Given a barcode-region sequence, outputs the sequence of the do-docamer barcode.
    This barcode (theoretically) consists of the concatentation of the two random hexamer sequences contained in the ligation oligo.
    However errors in sequences and ligation oligo production can mean that the random nucleotides are not always at the expected position.
    This function uses the known sequence of the spacers (which bound each of the two N6s to their 5') to deduce the random sequences.
    Returns a list of four numbmers, giving the start and stop positions of N1 and N2 respectively.
    """
    spcr = "GTCGTGAT"
    
    if "N" in bcseq and inputargs['allowNs'] == False:    # ambiguous base-call check 
      counts['getbarcode_fail_N'] += 1
      return
    
    # define expected positions
    first = bcseq[0:8]
    second = bcseq[14:22]   
    
    # if both spacers present and correct, return sequence from expected sites
    if first == spcr and second == spcr:
      counts['getbarcode_pass_exactmatch'] += 1
      return [8, 14, 22, 28]
    
    # otherwise look throughout the entire sequence for the presence of two spacers
    else:
      
      # pad the 5' of the sequence to allow for frame-shifts that result in incomplete primary spacer
      pad = 4
      bcseq = ("X"*pad) + bcseq 
      
      # search for all instances of the spacer, allowing for 2 substitutions first
        # failing that, search again, allowing for 2 substitutions OR (1 deletion or 1 insertion)
      err_spcrs = regex.findall("(GTCGTGAT){1s<=2}", bcseq) 
      positions = [bcseq.find(x) for x in err_spcrs]
      lens = [len(x) for x in err_spcrs]
      
      if len(positions) <> 2:
        err_spcrs = regex.findall("(GTCGTGAT){2i+2d+1s<=2}", bcseq) 
        positions = [bcseq.find(x) for x in err_spcrs]
        lens = [len(x) for x in err_spcrs]
        
        
        if len(positions) <> 2:
          counts['getbarcode_fail_not2spacersfound'] += 1
          return
      
      else:
        # if only two matches, first random seq runs from end of first spacer to start of next
          # second random seq just consists of the six bases following the second spacer
        n1pos = positions[0]+lens[0]
        n1len = positions[1] - n1pos
        n1end = n1pos + n1len
        
        n2pos = positions[1]+lens[1]
        n2len = 6
        n2end = n2pos + n2len

        if n1len <= 3:
          counts['getbarcode_fail_n1tooshort'] += 1
          return
        elif n1len >= 9:
          counts['getbarcode_fail_n1toolong'] += 1
          return
        elif n2end > len(bcseq):
          counts['getbarcode_fail_n2pastend'] += 1
          return          
        elif n1len == 6:
          counts['getbarcode_pass_fuzzymatch_rightlen'] += 1
        elif n1len in [4,5]:
          counts['getbarcode_pass_fuzzymatch_short'] += 1
        elif n1len >= 7:
          counts['getbarcode_pass_fuzzymatch_long'] += 1
        else:
          counts['getbarcode_pass_other'] += 1
          
        return [max(0,n1pos-pad), n1end-pad, n2pos-pad, n2end-pad]

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
    most_common_dcrseqs = [ds for ds, c in dcrseq_count.iteritems() if c == max_value]
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



def cluster_dcr_barcodes(counter, threshold):
    """
    Takes a counter {barcode:x, barcode:x, ...} of the number of copies of each barcode associated with a dcr
    Returns a counter {barcode_cluster:number_copies, ...} after combining barcodes that are below a threshold
    """

    barcodecluster_count = coll.Counter()
    if len(counter) == 1:
        # only one barcode is associated with this dcr, no clustering needed
        return counter

    track_barcodes_processed = coll.Counter()  # barcode:0/1 -> 1 if we have already incorporated the bc into a cluster
    for bc1, size1 in counter.most_common():
        if track_barcodes_processed[bc1] == 0:
            track_barcodes_processed[bc1] = 1
            barcodecluster_count[bc1] += size1
            for bc2, size2 in counter.most_common():
                if track_barcodes_processed[bc2] == 0:
                    if are_barcodes_equivalent(bc1, bc2, threshold):
                        track_barcodes_processed[bc2] = 1
                        barcodecluster_count[bc1] += size2

    return barcodecluster_count


def are_barcodes_equivalent(bc1, bc2, threshold):
    if lev.distance(bc1, bc2) <= threshold:
        return 1
    else:
        return 0


def makecounter():
    return coll.Counter()


def collapsinate(barcode_quality_parameters,
                 lev_threshold,
                 barcode_distance_threshold,
                 infile, 
                 outpath,
                 file_id):
  
    ###########################################
    ############# READING DATA IN #############
    ###########################################        
        
    # Check whether file appears to contain suitable verbose Decombinator output for collapsing
    if inputargs['dontcheckinput'] == False:
      if check_dcr_file(infile) != True:
        print "Please check that file contains suitable Decombinator output for collapsing."
        print "Alternatively, disable the input file sanity check by changing the \'dontcheckinput\' flag, i.e. \'-di True\'"
        sys.exit()
    
        
    print 'Reading data in...'
    t0 = time()
    barcode_dcretc = coll.defaultdict(list)
    input_dcr_counts = coll.Counter()
    inhandle = opener(infile, 'r')
    
    for line in inhandle:
        counts['readdata_input_dcrs'] += 1
        fields = line.rstrip('\n').split(', ')

        bc_locs = get_barcode(fields[8])        # barcode locations
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
            barcode_dcretc[barcode].append(dcretc)
        else:
            counts['readdata_fail_no_bclocs'] += 1

    inhandle.close()
    
    counts['readdata_barcode_dcretc_keys'] = len(barcode_dcretc.keys())
    counts['number_input_unique_dcrs'] = len(input_dcr_counts.keys())
    counts['number_input_total_dcrs'] = sum(input_dcr_counts.values())
      
    t1 = time()
    print '  ', round(t1-t0, 2), 'seconds'
    counts['time_readdata_s'] = t1
    
    ##########################################
    ############# TCR COLLAPSING #############
    ##########################################        
    
    print 'Collapsing by barcode...'
    t0 = time()
    
    dcr_barcodecounter = coll.defaultdict(makecounter)
    
    for barcode, dcretc_list in barcode_dcretc.iteritems():

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
    print '  ', round(t1-t0, 2), 'seconds'
    counts['time_tcrcollapsing_s'] = t1
    

    ##############################################
    ############# BARCODE CLUSTERING #############
    ##############################################          

    print 'Clustering barcodes...'
    t0 = time()
    
    dcr_originalcopies = coll.Counter()         # Stores the final collapsed DCRs with estimated frequencies (i.e. # clustered barcodes)
    barcode_duplication = coll.Counter()        # Stores the frequency with which each barcode appeared in input file (number lines)
    
    for dcr, barcode_count in dcr_barcodecounter.iteritems():

        clustered_bc_dict = cluster_dcr_barcodes(barcode_count, barcode_distance_threshold)
        dcr_originalcopies[dcr] = len(clustered_bc_dict)
        for bc, c in clustered_bc_dict.iteritems():
            barcode_duplication[bc] += c
    
    counts['number_output_unique_dcrs'] = len(dcr_originalcopies)
    counts['number_output_total_dcrs'] = sum(dcr_originalcopies.values())
    
    t1 = time()
    counts['time_bcclustering_s'] = t1
    print '  ', round(t1-t0, 2), 'seconds'

    outfile = outpath + file_id + suffix
    outhandle = open(outfile, 'w')
    print 'Writing to output file', outfile, '...'
    
    for dcr, copies in dcr_originalcopies.iteritems():
        print >> outhandle, ', '.join([dcr, str(copies)])
    outhandle.close()
    
    if inputargs['dontgzip'] == False:  # Gzip output file
        print "Compressing output to", outfile + ".gz ..."
      
        with open(outfile) as inf, gzip.open(outfile + '.gz', 'wb') as outf:
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
        for bc, copies in barcode_duplication.iteritems():
            print >> outhandle, ','.join([bc, str(copies)])
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
      
      summstr = "Property,Value\nDirectory," + os.getcwd() + "\nInputFile," + inputargs['infile'] + "\nOutputFile," + counts['outfilenam'] \
        + "\nDateFinished," + date + "\nTimeFinished," + strftime("%H:%M:%S") + "\nTimeTaken(Seconds)," + str(round(counts['time_taken_total_s'],2)) + "\n\n"
      
      for s in ['extension', 'dontgzip', 'allowNs', 'dontcheckinput', 'barcodeduplication', 'minbcQ', 'bcQbelowmin', 'bcthreshold', \
        'lenthreshold', 'percentlevdist', 'avgQthreshold']:
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

      print >> summaryfile, summstr 
      summaryfile.close()

