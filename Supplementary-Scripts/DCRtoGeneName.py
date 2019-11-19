# DCRtoGeneName.py 
# March 2016, Jamie Heather, UCL

# Take a file containing a DCR identifier in its first 5 comma-delimited fields (e.g. n12, freq or dcrcdr3) and convert the first two to the actual gene name

import argparse
import re
import sys
import collections as coll
import os
import urllib.request
import gzip

# Fix - add an optional outputting of functionality? Need to use .translate files instead of tag files.

def args():
    """args(): Obtains command line arguments which dictate the script's behaviour"""

    # Help flag
    parser = argparse.ArgumentParser(
        description='Convert Decombinator V and J indexes to their IMGT equivalent gene names.')
    # Add arguments
    parser.add_argument(
        '-c', '--chain', type=str, help='TCR chain (a/b/g/d)', required=False)
    parser.add_argument(
        '-in', '--infile', type=str, help='Optional input file', required=False)
    parser.add_argument(
        '-dz', '--dontgzip', type=bool, help='Stop the output files automatically being compressed with gzip (True/False)', required=False, default=False)
    parser.add_argument(
        '-ex', '--extension', type=str, help='Specify the file extension of the output DCR file. Default is to append \"tr_\" to existing.', required=False)
    parser.add_argument(
        '-pf', '--prefix', type=str, help='Specify the prefix of the output file. None set by default', required=False, default="")
    parser.add_argument(
        '-tg', '--tags', type=str, help='Specify which Decombinator tag set to use (extended or original). Default = extended', required=False, default="extended")
    parser.add_argument(
        '-sp', '--species', type=str, help='Specify which species TCR repertoire the data consists of (human or mouse). Default = human', required=False, default="human")
    parser.add_argument(
        '-v', '--v_genes', type=int, help='Numeric Decombinator index of V gene', required=False, nargs='+')    
    parser.add_argument(
        '-j', '--j_genes', type=int, help='Numeric Decombinator index of J gene', required=False, nargs='+')        
    parser.add_argument(
        '-tfdir', '--tagfastadir', type=str, help='Path to folder containing TCR FASTA and Decombinator tag files, for offline analysis. \
        Default = \"Decombinator-Tags-FASTAs\".', required=False, default="Decombinator-Tags-FASTAs")

    return parser.parse_args()

def commas(line):
    """ Return the location of all the commas in an input string or line """
    return [m.start() for m in re.finditer(',', line)] 

def read_tcr_file(species, tagset, gene, filetype, expected_dir_name):
  """
  Reads in the associated data for the appropriate TCR locus from the ancillary files (hosted in own repo)
  :param species: human or mouse
  :param tagset: original or extended
  :param gene: V or J
  :param filetype: tag/fasta/translate/cdrs
  :param expected_dir_name: (by default) Decombinator-Tags-FASTAs
  :return: the opened file (either locally or remotely)
  """
  # Define expected file name
  expected_file = species + "_" + tagset + "_" + "TR" + chain.upper() + gene.upper() + "." + filetype

  # First check whether the files are available locally (in pwd or in bundled directory)
  if os.path.isfile(expected_file):
    fl = expected_file

  elif os.path.isfile(expected_dir_name + os.sep + expected_file):
    fl = expected_dir_name + os.sep + expected_file

  else:
    try:
      fl = "https://raw.githubusercontent.com/innate2adaptive/Decombinator-Tags-FASTAs/master/" + expected_file
      urllib.request.urlopen(fl)  # Request URL, see whether is found
      fl = urllib.request.urlretrieve(fl)[0]

    except Exception:
      print("Cannot find following file locally or online:", expected_file)
      print("Please either run Decombinator with internet access, or point Decombinator to local copies of the tag and FASTA files with the \'-tf\' flag.")
      sys.exit()

    # Return opened file, for either FASTA or tag file parsing
  return fl

def import_tcr_info(inputargs):
    """ import_tcr_info: Gathers the required TCR chain information for Decombining """
      
    # Get chain information
    global chainnams, chain, counts
    counts = coll.Counter()
    chainnams = {"a": "alpha", "b": "beta", "g": "gamma", "d": "delta"}
      
    nochain_error = "TCR chain not recognised. \n \
    Please either include (one) chain name in the file name (i.e. alpha/beta/gamma/delta),\n \
    or use the \'-c\' flag with an explicit chain option (a/b/g/d, case-insensitive)."
    
    if inputargs['chain']:
      if inputargs['chain'].upper() in ['A', 'ALPHA', 'TRA', 'TCRA']:
        chain = "a" 
      elif inputargs['chain'].upper() in ['B', 'BETA', 'TRB', 'TCRB']:
        chain = "b" 
      elif inputargs['chain'].upper() in ['G', 'GAMMA', 'TRG', 'TCRG']:
        chain = "g" 
      elif inputargs['chain'].upper() in ['D', 'DELTA', 'TRD', 'TCRD']:
        chain = "d" 
      else:
        print(nochain_error)
        sys.exit()
    else:
      # If no chain provided, try and infer from filename
      inner_filename_chains = [x for x in chainnams.values() if x in inputargs['infile'].lower()]
      if len(inner_filename_chains) == 1:
        chain = inner_filename_chains[0][0]  
      else:
        print(nochain_error)
        sys.exit()
      
    print('Importing TCR', chainnams[chain], 'gene sequences...')

    # First check that valid tag/species combinations have been used
    if inputargs['tags'] == "extended" and inputargs['species'] == "mouse":
      print("Please note that there is currently no extended tag set for mouse TCR genes.\n \
      Script will now switch the tag set in use from \'extended\' to \'original\'.\n \
      In future, consider editing the script to change the default, or use the appropriate flags (-sp mouse -tg original).")
      inputargs['tags'] = "original"
    
    if inputargs['tags'] == "extended" and ( chain == 'g' or chain == 'd' ):
      print("Please note that there is currently no extended tag set for gamma/delta TCR genes.\n \
      Script will now switch the tag set in use from \'extended\' to \'original\'.\n \
      In future, consider editing the script to change the default, or use the appropriate flags.")
      inputargs['tags'] = "original"

    # FIX check tag set.

    # Check species information
    if inputargs['species'] not in ["human", "mouse"]:
      print("Species not recognised. Please select either \'human\' (default) or \'mouse\'.\n \
      If mouse is required by default, consider changing the default value in the script.")
      sys.exit()    
      
    # Look for tag and V/J fasta and tag files: if these cannot be found in the working directory, source them from GitHub repositories
      # Note that fasta/tag files fit the pattern "species_tagset_gene.[fasta/tags]"
      # I.e. "[human/mouse]_[extended/original]_TR[A/B/G/D][V/J].[fasta/tags]"
    
    for gene in ['v', 'j']:
      # Get names from tag files
      tag_file = read_tcr_file(inputargs['species'], inputargs['tags'], gene, "tags", inputargs['tagfastadir'])  # get tag data
      tag_data = open(tag_file, "rt")     
      globals()[gene+"_nams"] =  get_gene_names(tag_data)
      tag_data.close()
    
def get_gene_names(vj_file):
    """Read V names in from tag file"""
    nams = [] 
    for line in vj_file:
      elements = line.rstrip("\n").split('|')
      nams.append(elements[1].split('*')[0]) # Adds elements in first column iteratively
    return nams

def sort_permissions(fl):
    # Need to ensure proper file permissions on output data
      # If users are running pipeline through Docker might otherwise require root access
    if oct(os.stat(fl).st_mode)[4:] != '666':
      os.chmod(fl, 0o666)

if __name__ == '__main__':

    inputargs = vars(args())
    counts = coll.Counter()

    chainnams = {"a": "alpha", "b": "beta", "g": "gamma", "d": "delta"}

    import_tcr_info(inputargs)
  
    if inputargs['v_genes'] == None and inputargs['j_genes'] == None and inputargs['infile'] == None:
      print("Script requires something to convert: either a dcr/n12/freq/dcrcdr3 file (use the -in flag) or a V and/or J gene (use the -v or -j flags).")
      sys.exit() 
      
    # If users specify V or J genes in the command line, then print these to stdout
    if inputargs['v_genes'] != None:
      if max(inputargs['v_genes']) > len(v_nams):
        print("Input V index detected that exceeds the range in use for this tag set (" \
        + str(max(inputargs['v_genes'])) + " vs " + str(len(v_nams)) + ")")
        sys.exit()
      print("\nInput V indexes:\t" + str(inputargs['v_genes']))
      print("Converted gene names:\t" + str([v_nams[x] for x in inputargs['v_genes']]))
    if inputargs['j_genes'] != None:
      if max(inputargs['j_genes']) > len(j_nams):
        print("Input J index detected that exceeds the range in use for this tag set (" \
        + str(max(inputargs['j_genes'])) + " vs " + str(len(j_nams)) + ")")
        sys.exit()
      print("\nInput J indexes:\t" + str(inputargs['j_genes']))
      print("Converted gene names:\t" + str([j_nams[x] for x in inputargs['j_genes']]))

    if inputargs['infile'] != None:
      # Set file openers
      if inputargs['infile'].endswith('.gz'):
        in_opener = gzip.open
      else:
        in_opener = open
      
      samplenam = inputargs['infile'].split('.')[0]
      if inputargs['extension'] != None:
        new_extension = inputargs['extension']
      else:
        previous_extension = inputargs['infile'].split('.')[1]
        new_extension = "tr_" + previous_extension
      
      if inputargs['dontgzip'] == True:
        out_opener = open
        
        outfilename = inputargs['prefix'] + samplenam + '.' + new_extension
        if outfilename.endswith('.gz'):
          outfilename = outfilename[:-3]
      else:
        out_opener = gzip.open
        outfilename = inputargs['prefix'] + samplenam + '.' + new_extension + '.gz'
      
      with in_opener(inputargs['infile'], 'rt') as infile, out_opener(outfilename, 'wt') as outfile:
        
        print("Reading from", inputargs['infile'], "and writing to", outfilename)
        for line in infile:
          counts['in_line'] += 1
          comma =  commas(line)
          if counts['in_line'] < 5:
            if len(comma) >= 4:
              try:
                int(line[:comma[0]])
              except:
                print("File does not appear to contain Decombinator index data.")
                sys.exit()
          
          v_num = int(line[:comma[0]])
          j_num = int(line[comma[0]+1:comma[1]])

          if v_num not in range(len(v_nams)) or j_num not in range(len(j_nams)):
            print("Error detected in the following line, as gene index does not fall in expected range. Ensure correct species, chain and tag flags are set.")
            print(line)
            sys.exit()
          else:
            v = v_nams[v_num]
            j = j_nams[j_num]
            
            outline = v + ", " + j + line[comma[1]:]
            # from IPython import embed
            # embed()
            outfile.write(outline)

      sort_permissions(outfilename)
    




