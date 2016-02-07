# ligTCRtranslateCDR3.py 
# James M. Heather, February 2016, UCL

##################
### BACKGROUND ###
##################

# Built on CDR3ulator.py v2

# Take any decombined data and output the functional CDR3s only
  # CDR3s must be: in-frame; lacking-stop codons, and run from a conserved cysteine to FGXG motif (or appropriate alternatives)
# Originally built on PlotDCR.py (itself modified from dcr v1.4), then algorithm swapped
  # Now uses CDR3 detection based on Katharine's functions.py script
  # Uses defined conserved cysteine residues, and allows for atypical TRAJ CDR3-ending motifs
  
# Makes use of functions originally developed in Katharine Best's functions.py script and Niclas Thomas' Decombinator (v1.2)
# Note that this version provides the capabilities to generate CDR3s from all genes covered in the extended Decombinator tags
  # Many of these genes will never generate CDR3s from this code regardless of whether they're included in CDR3s
  # This is typically due to V genes that lack the conserved C residue, or contain stop codons upstream of it.
  # These have been left in, both to provide a single location that describes the whole heterogeneity of the prototypical alpha/beta genes

# New in v2:
  # Have an option to output a file of non-functional rearrangements 
  # Provides an option to turn off statistics standard out results
    # This now includes the percentages of the different reasons for being assigned non-functional

##################
###### INPUT #####
##################

# Takes any text file in comma-space (", ") delimited decombinator format
  # 5-part TCR identifier must come first, followed by any additional fields, such as frequency

# BIG FIX - sort documentation

# BIG FIX - make translation files for original/mouse/gammadelta


# Note that in addition to the same FASTA files that Decombinator makes use of, this script requires additional '.translate' files
  # These contain four comma-delimited fields, which allow for the correct translation of TCR sequences from DCR indexes
  # Those fields are: Gene name, conserved position (of cysteine or FGXG motif), conserved residue (when the C/FGXG differs), and designated functionality

##################
##### OUTPUT #####  
##################

# Users have options of two output, comma-delimited formats:
  # '.cdr3' files, which consist of the unique productively-rearranged CDR3s from the original file, with frequencies
  # '.dcrcdr3' files, which contains the five-part Decombinator index before the CDR3 it encodes and its frequency
  # The choice of which file format is used is decided by altering the 'inputargs['dcroutput']' variable
  # Note that output files will likely be shorter than the input file, due to multiple TCR sequences encoding the same CDR3
# Users also have the option to include the final three residues of the CDR3 in the output, or to stop at the phenylalanine 

# The script also outputs some simple statistics about the number of productive rearrangements, and the frequency of different functionality genes
# NB: See IMGT for further information on definitions of functionality and productivity:
  # http://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html
  # Essentially, 'functionality' refers to the predicted ability of certain germline genes to contribute to working receptors
  # 'Productive' refers to the predicted ability of a given rearrangement to encode a working receptor chain
  
##################
#### PACKAGES ####  
##################

from __future__ import division
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
import argparse
import string
import re
import sys
import collections as coll
import os
import urllib2
import warnings

# Supress Biopython translation warning when translating sequences where length % 3 != 0
warnings.filterwarnings("ignore") 
   
   # BIG FIXES - HUMAN/MOUSE, A/B/G/D, NEW SOURCE TAGS ETC
   
   
###################
#### FUNCTIONS ####
###################

def args():
    """args(): Obtains command line arguments which dictate the script's behaviour"""

    # Help flag
    parser = argparse.ArgumentParser(
        description='Translate and extract CDR3 sequences from Decombinator classifier files.')
    # Add arguments
    parser.add_argument(
        '-in', '--infile', type=str, help='File containing 5 part Decombinator indexes, (with/without frequencies)', required=True)
    parser.add_argument(
        '-c', '--chain', type=str, help='TCR chain (a/b/g/d)', required=False)
    parser.add_argument(
        '-sp', '--species', type=str, help='Specify which species TCR repertoire the data consists of (human or mouse). Default = human', required=False, default="human")
    parser.add_argument(
        '-tg', '--tags', type=str, help='Specify which Decombinator tag set to use (extended or original). Default = extended', required=False, default="extended")
    parser.add_argument(
        '-s', '--suppresssummary', type=bool, help='Output summary data (True/False)', required=False, default=False)
    parser.add_argument(
        '-dz', '--dontgzip', type=bool, help='Stop the output FASTQ files automatically being compressed with gzip (True/False)', required=False, default=False)
    parser.add_argument(
        '-dc', '--dontcount', type=bool, help='Show the count (True/False)', required=False, default=False)
    parser.add_argument(
        '-ex', '--extension', type=str, help='Specify the file extension of the output DCR file. Default = \"n12\"', required=False, default="n12")
    parser.add_argument(
        '-do', '--dcroutput', type=bool, help='Optionally include Decombinator TCR index along with the CDR3 sequence and frequency. Default = False', \
        required=False, default=False)
    parser.add_argument(
        '-tfdir', '--tagfastadir', type=str, help='Path to folder containing TCR FASTA and Decombinator tag files, for offline analysis. \
        Default = \"Decombinator-Tags-FASTAs\".', required=False, default="Decombinator-Tags-FASTAs")    
    parser.add_argument(
        '-gxg', '--includeGXG', type=bool, help='Optionally include the \"GXG\" motif following the conserved phenylalanine residue that terminates the CDR3 region. Defaut = False', required=False, default=False)
    parser.add_argument(
        '-np', '--nonproductive', type=bool, help='Optionally output an additional file containing the non-productive TCR rearrangements. Default =  False', required=False, default=False)
    
    # fix - add option to change suffixes (productive, non-prod and dcr containing)
    
    return parser.parse_args()

def findfile(testfile):
    """ Check whether file is present at given path """
    try:
        testopen = open(str(filename),"rU")
        testopen.close()
    except:
        print 'Cannot find the specified input file. Please try again'
        sys.exit()


def read_tcr_file(species, tagset, gene, filetype, expected_dir_name):
  """ Reads in the FASTA and tag data for the appropriate TCR locus """
  
  # Define expected file name
  expected_file = species + "_" + tagset + "_" + "TR" + chain.upper() + gene.upper() + "." + filetype

  # First check whether the files are available locally (in pwd or in bundled directory)
  if os.path.isfile(expected_file):
    fl = expected_file
    fl_opener = open
  elif os.path.isfile(expected_dir_name + os.sep + expected_file):
    fl = expected_dir_name + os.sep + expected_file
    fl_opener = open
  else:
    try:
      fl = "https://raw.githubusercontent.com/JamieHeather/Decombinator-Tags-FASTAs/master/" + expected_file
      urllib2.urlopen(urllib2.Request(fl))      # Request URL, see whether is found
      fl_opener = urllib2.urlopen
    except:
      print "Cannot find following file locally or online:", expected_file
      print "Please either run Decombinator with internet access, or point Decombinator to local copies of the tag and FASTA files with the \'-tf\' flag."
      sys.exit()
  
  # Return opened file, for either FASTA or tag file parsing
  return fl_opener(fl)

def import_gene_information(inputargs):            
    """ Obtains gene-specific information for translation """
    
    # Runs first: reads in V and J gene sequence and name data (from fasta files)
      # and positions of conserved cysteine residues in V genes (from separate files)
    
    # If files cannot be found in local directory, script looks for them online at GitHub
    
    # NB that a number of psuedogenes have no officially designated conserved C (or indeed a 3' C at all)
      # Where possible, the nearest suitable C residue is used, where not an arbitrary position of 0 is given
      # Somewhat moot, as most psuedogenes contain a number of stop codons and thus cannot produce productive rearrangements 
      
    # First check that valid tag/species combinations have been used
    if inputargs['tags'] == "extended" and inputargs['species'] == "mouse":
      print "Please note that there is currently no extended tag set for mouse TCR genes.\n \
      Decombinator will now switch the tag set in use from \'extended\' to \'original\'.\n \
      In future, consider editing the script to change the default, or use the appropriate flags (-sp mouse -tg original)."
      inputargs['tags'] = "original"
    
    if inputargs['tags'] == "extended" and ( chain == 'g' or chain == 'd' ):
      print "Please note that there is currently no extended tag set for gamma/delta TCR genes.\n \
      Decombinator will now switch the tag set in use from \'extended\' to \'original\'.\n \
      In future, consider editing the script to change the default, or use the appropriate flags."
      inputargs['tags'] = "original"

    # Check species information
    if inputargs['species'] not in ["human", "mouse"]:
      print "Species not recognised. Please select either \'human\' (default) or \'mouse\'.\n \
      If mouse is required by default, consider changing the default value in the script."
      sys.exit()    
      
    # Look for tag and V/J fasta and cysteine position files: if these cannot be found in the working directory, source them from GitHub repositories
      # Note that fasta/tag files fit the pattern "species_tagset_gene.[fasta/tags]"
      # I.e. "[human/mouse]_[extended/original]_TR[A/B/G/D][V/J].[fasta/tags]"
    
    for gene in ['v', 'j']:
      # Get FASTA data
      fasta_file = read_tcr_file(inputargs['species'], inputargs['tags'], gene, "fasta", inputargs['tagfastadir'])  
      globals()[gene+"_genes"] = list(SeqIO.parse(fasta_file, "fasta"))
      fasta_file.close()
      globals()[gene+"_regions"] = [str(string.upper(item.seq)) for item in globals()[gene+"_genes"]]  
      globals()[gene+"_names"] = [str(string.upper(item.id).split("|")[1]) for item in globals()[gene+"_genes"]]  
      
      # Get conserved translation residue sites and functionality data
      translation_file = read_tcr_file(inputargs['species'], inputargs['tags'], gene, "translate", inputargs['tagfastadir'])  
      translate_data = [x.rstrip() for x in list(translation_file)]
      translation_file.close()
      globals()[gene+"_translate_position"] = [int(x.split(",")[1]) for x in translate_data]
      globals()[gene+"_translate_residue"] = [x.split(",")[2] for x in translate_data]
      globals()[gene+"_functionality"] = [x.split(",")[3] for x in translate_data]
      
    return v_regions, j_regions, v_names, j_names, v_translate_position, v_translate_residue, j_translate_position, j_translate_residue,\
      v_functionality, j_functionality

def get_cdr3(dcr, chain, vregions, jregions, vtranslate_pos, vtranslate_res, jtranslate_pos, jtranslate_res):
    """ Checks the productivity of a given DCR-assigned rearrangement 
    Returns a 1 if productive, 0 if not """
    
    # NB: A productively rearranged receptor does not necessarily mean that it is the working receptor used in a cell
      # It could be a silenced chain that isn't used, or could have inactivating mutations upstream of the sequenced region
    
    # 0.5 Set up check variables
    
    # Boolean productivity checks that CDR3s must pass 
    in_frame = 0
    no_stop = 0
    found_c = 0
    found_fgxg = 0
    
    # CDR3-defining positions 
    start_cdr3 = 0
    end_cdr3 = 0 
    
    # 1. Rebuild whole nucleotide sequence from Decombinator assignment
    classifier_elements = dcr.split(', ')
    v = int(classifier_elements[0])
    j = int(classifier_elements[1])
    vdel = int(classifier_elements[2])
    jdel = int(classifier_elements[3])
    ins_nt = classifier_elements[4]

    if vdel == 0:
      v_used = vregions[v]
    else:
      v_used = vregions[v][:-vdel]
      
    j_used = jregions[j][jdel:]

    nt = ''.join([v_used, ins_nt, j_used])
    
    # 2. Translate
    aa = str(Seq(nt, generic_dna).translate())

    # 3. Check whether whole rearrangement is in frame
    if (len(nt)-1) % 3 == 0:
      in_frame = 1
    else:
      if '*' in aa:
        return "OOF_with_stop"
      else:
        return "OOF_without_stop"

    # 4. Check for stop codons in the in-frame rearrangements
    if '*' not in aa:
      no_stop = 1
    else:
      return "IF_with_stop"

    # 5. Check for conserved cysteine in the V gene
    cdr3_start = vtranslate_pos
    cdr3_c = vtranslate_res
    
    if aa[vtranslate_pos[v]-1] == vtranslate_res[v]:
      found_c = 1
      start_cdr3 = vtranslate_pos[v]-1
    else:
      return "No_conserved_cysteine"
    
    # 5.5 Having found conserved cysteine, only need look downstream to find other end of CDR3
    downstream_c = aa[start_cdr3:]

    # 6. Check for presence of FGXG motif (or equivalent)
    site = downstream_c[jtranslate_pos[j]:jtranslate_pos[j]+4]
    
    if re.findall(jtranslate_res[j], site):
      if inputargs['includeGXG'] == True: 
        end_cdr3 = len(downstream_c) - 11 + start_cdr3 + 4 
      else:
        end_cdr3 = len(downstream_c) - 11 + start_cdr3 + 1
    
    else:
      return "No_conserved_FGXG"
    
    return aa[start_cdr3:end_cdr3]    

####################################
######## CHECK INPUT FILES #########
####################################

if __name__ == '__main__':
    
    # Get parameters 
    inputargs = vars(args())
    counts = coll.Counter()
    
    if inputargs['infile'].endswith('.gz'):
        opener = gzip.open
    else:
        opener = open    
      
    #  Get chain information
    if not inputargs['chain']:
        # If chain not given, try and infer from input file name 
        chaincheck = [x for x in ["alpha", "beta", "gamma", "delta"] if x in inputargs['infile'].lower()]
        if len(chaincheck) == 1:
            chain = chaincheck[0]
        else:
            print "TCR chain not recognised. Please choose from a/b/g/d (case-insensitive)."
            sys.exit()          
    else:
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
    
    suffix = "." + inputargs['extension']

    filename = inputargs['infile']
    findfile(filename)
    
    ####################################
    ########## EXTRACT CDR3s ###########
    ####################################

    v_regions, j_regions, v_names, j_names, v_translate_position, v_translate_residue, j_translate_position, j_translate_residue, \
      v_functionality, j_functionality = import_gene_information(inputargs)

    infile = open(filename, "rU")

    counts['line_count'] = 0

    # Count non-productive rearrangments
    fail_count = coll.Counter()
    fails = ["OOF_with_stop", "OOF_without_stop", "IF_with_stop", "No_conserved_cysteine", "No_conserved_FGXG"]

    # Store and count CDR3s
    if inputargs['dcroutput'] == False:
      cdr3_count = coll.Counter()
    elif inputargs['dcroutput'] == True:
      dcr_cdr3_count = coll.Counter()

    np_cdr3_count = coll.Counter()

    ## Count the number of F, ORF and P genes used for both productive and non-productive rearrangements # FIX rm - change so adds to counts
    #pVf = 0
    #pVorf = 0
    #pVp = 0
    #pJf = 0
    #pJorf = 0
    #pJp = 0

    #npVf = 0
    #npVorf = 0
    #npVp = 0
    #npJf = 0
    #npJorf = 0
    #npJp = 0

    for line in infile:
      
      counts['line_count'] += 1
      
      comma = [m.start() for m in re.finditer(',', line)]           

      in_dcr = str(line[:comma[4]])
      
      cdr3 = get_cdr3(in_dcr, chain, v_regions, j_regions, v_translate_position, v_translate_residue, j_translate_position, j_translate_residue)
      
      frequency = int(line[comma[4]+2:].rstrip())
      
      dcr_cdr3 = in_dcr + ":" + cdr3
      
      v = int(line[:comma[0]])
      j = int(line[comma[0]+2:comma[1]])
      
      
          #info = import_gene_information(inputargs) # FIX - rm these old two lines if not used
      # Returns: 0 v_regions, 1 j_regions, 2 v_names, 3 j_names, 4 v_translate, 5 j_translate
    #v_regions, j_regions, v_names, j_names, v_translate_position, v_translate_residue, j_translate_position, j_translate_residue \
      #= import_gene_information(inputargs)

      
      
      
      if cdr3 not in fails:
        counts['prod_recomb'] += 1    
        productivity = "P"
        
        if inputargs['dcroutput'] == False:
          cdr3_count[cdr3] += frequency
        elif inputargs['dcroutput'] == True:
          dcr_cdr3_count[dcr_cdr3] += frequency
            
      else:  
        np_cdr3_count[dcr_cdr3] += frequency
        fail_count[cdr3] += 1
        productivity = "NP"
      
      # Count the number of number of each type of gene functionality (by IMGT definitions, based on prototypic gene) #fix RM
      counts[productivity + "_" + "V-" + v_functionality[v]] += 1
      counts[productivity + "_" + "J-" + j_functionality[j]] += 1
      
    if inputargs['dcroutput'] == True:

      outfilename = filename.split(".")[0]+".dcrcdr3"
      outfile = open(outfilename, "w")
        
      for x in dcr_cdr3_count:
        outtext = x + ", " + str(dcr_cdr3_count[x])
        print >> outfile, outtext
          
      infile.close()
      outfile.close()
      
    elif inputargs['dcroutput'] == False:

      outfilename = filename.split(".")[0]+".cdr3"
      outfile = open(outfilename, "w")

      for x in cdr3_count:
        outtext = x + ", " + str(cdr3_count[x])
        print >> outfile, outtext
          
      infile.close()
      outfile.close()

    if inputargs['nonproductive'] == True:  
      npfilename = filename.split(".")[0]+".np"
      npfile = open(npfilename, "w")

      for x in np_cdr3_count:
        outtext = x + ", " + str(np_cdr3_count[x])
        print >> npfile, outtext
      npfile.close()    
        
    NP_count = sum(fail_count.values())
            
    ###################
    ##### RESULTS #####
    ###################

    sys.exit()
    # fix output stats to summary file
  
    print "Reading", str(counts['line_count']), "Decombinator-assigned rearrangements from", str(filename) + ", and writing out to", str(outfilename), "\n"

    print '{0:,}'.format(counts['prod_recomb']), "productive rearrangements detected" 
    print "\tV gene usage:\t" + '{0:,}'.format(pVf), "F;\t" + '{0:,}'.format(pVorf), "ORF;\t" +  '{0:,}'.format(pVp), "P"
    print "\tJ gene usage:\t" + '{0:,}'.format(pJf), "F;\t" + '{0:,}'.format(pJorf), "ORF;\t" + '{0:,}'.format(pJp), "P\n"

    print '{0:,}'.format(NP_count), "non-productive rearrangements detected"
    print "\tV gene usage:\t" + '{0:,}'.format(npVf), "F;\t" + '{0:,}'.format(npVorf), "ORF;\t" +  '{0:,}'.format(npVp), "P"
    print "\tJ gene usage:\t" + '{0:,}'.format(npJf), "F;\t" + '{0:,}'.format(npJorf), "ORF;\t" + '{0:,}'.format(npJp), "P\n"

    print 'Non-productive rearrangement statistics:'
    print '\tOut of frame, no stop codon:\t {:.2%}'.format(fail_count['OOF_without_stop']/NP_count) + "\t(" + str(fail_count['OOF_without_stop']) + ")"
    print '\tOut of frame, with stop codon:\t {:.2%}'.format(fail_count['OOF_with_stop']/NP_count) + "\t(" + str(fail_count['OOF_with_stop']) + ")"
    print '\tIn frame, with stop codon:\t {:.2%}'.format(fail_count['IF_with_stop']/NP_count) + "\t(" + str(fail_count['IF_with_stop']) + ")"
    print '\tNo conserved cysteine at start:\t {:.2%}'.format(fail_count['No_conserved_cysteine']/NP_count) + "\t(" + str(fail_count['No_conserved_cysteine']) + ")"
    print '\tNo conserved FGXG at end:\t {:.2%}'.format(fail_count['No_conserved_FGXG']/NP_count) + "\t(" + str(fail_count['No_conserved_FGXG']) + ")"

    if (counts['prod_recomb'] + NP_count) <> line_count:
      print "\nError detected: Sum of productive and non-productive sorted sequences does not equal total number of input sequences"
    