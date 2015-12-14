# vDCR v1.5 : "verbose Decombinator"
# James M. Heather, August 2014, UCL

##################
### BACKGROUND ###
##################

# This script reads in T-cell receptor (TCR) deep-sequencing data and assigns them a Decombinator (DCR) index
# This script is itself a modified version of Decombinator v1.4 (by Dr. Niclas Thomas)
# However I have modified it to allow downstream error-correction by making use of random DNA barcodes incorporated before PCR
# This version additionally incorporates other changes, including increased stringency filters and different input methods
# This version also provides the option to use a new set of Decombinator V/J gene tags (see 'tags' variable)
  # These new extended tag sets are able to identity all IMGT-proscribed prototypic genes, regardless of functional status
# For details on the official version of Decombinator, see:
  # http://dx.doi.org/10.1093/bioinformatics/btt004 (Thomas et al 2013 publication)
  # https://github.com/uclinfectionimmunity/Decombinator/ (GitHub repository)
# NB Many of the changes made in this version have been included in the official v2 Decombinator releases

##################
###### INPUT #####
##################

# Takes fastq files in which the first 12 nucleotides are made up of random barcode sequence introduced before PCR amplification
# Such fastq files can be produced by DualIndexDemultiplexing.py or AddAllR2Hex.py
# Only works for human alpha and beta chain TCR sequence data

# Takes command line input, have to specify which chain to Decombine (a/b), e.g.:
  # python vDCR.py FILE1.fq a

# Users can also pipe (likely gzipped) data in, e.g.:
  # zcat FILE2.fq.gz | python vDCR.py FILE2.fq.gz b
  # NB filename is required again after pipe in order to correctly name output files

# Note that this version also doesn't necessarily require the tag and V/J fasta files to be in the same directory as the script
  # If it cannot find them locally, the script will instead attempt to download them from GitHub

##################
##### OUTPUT #####  
##################

# Produces a '.n12' file, which is a standard comma-delimited Decombinator output file with several additional fields:
  # V index, J index, # V deletions, # J deletions, insert, ID, tcr sequence, tcr quality, barcode sequence, barcode quality
# Also outputs some statistics on to stdout
# .n12 files can be processed by CollapseTCRs.py, to mitigate for PCR and sequencing errors and amplification

# NB By default the omitN option (which skips reads that contain an ambiguous base call in any position) is set to False
  # Instead I recommend leaving these reads, and then only filtering out those that contain Ns in barcode or inter-tag region
  # e.g.: for i in *n12; do awk '{FS=", "} {if ($7~/N/ || $9~/N/){next} else{ if (length($7) < 130){print}}}' $i > temp && mv temp $i; done
  # This one-liner also removes TCRs that contain erroneous assignations with over-long inter-tag regions

# Output files are named based on input filename, using '.' dots to remove suffixes, thus these must not be included in filenames
# Performs no CDR3 translation or plotting; these functionalities are available in other scripts in the suite


##################
#### ANALYSIS ####
##################


import traceback    
import sys          
import os
import platform
import urllib2

tags = "extended"
  # Users can specify here which set of Decombinator tags they wish to use 
  # Options are "original" and "extended"
  # "original" tags:
    # Those that were included in original Thomas et al 2013 paper
    # Tags generated for the prototypic allele of each gene designated functional ('F') by IMGT 
  # "extended" tags:
    # Includes tags for the prototypic allele of every gene, including 'P' and 'ORF' IMGT-designations
    # This allows one to search for more TCR recombinations that exist but aren't necessarily functional
  # NB Output files will be prefixed differently according to which tag set is used:
    # Files decombined with original tags are prefixed 'vDCR_"
    # Files decombined with extended tags are prefixed 'vDCRe_"
  # Note that extended have had the new genes added after those that were in the original set
    # I.e. Each particular gene index from original tags will be the same in the new tags
    # However as some of the actual tag sequences have changed, the different tag sets might not find the exact same TCRs
    

def analysis( Sequence_Reads, results, chain, with_statistics=True, with_reverse_complement_search=True, omitN=False):
    import numpy as np
    import decimal as dec
    import string
    import operator as op
    import collections as coll
    from Bio import SeqIO
    from acora import AcoraBuilder
    from time import time, clock
    from string import Template
    from operator import itemgetter, attrgetter
    import Levenshtein as lev

    imposs_deletion = 0
    tag_overlap = 0
    
    results_file = open( str(results)+'.n12', "w")      

    Nseqs = 0

    if chain=="alpha":

        # Do not change - V tags are split at position 10, J at position 10, to look for half tags if no full tag is found.
        v_half_split, j_half_split = [10,10] 

        ################

        print 'Commencing analysis on a total of', len(Sequence_Reads), 'file(s)'

        ################
        print ('Importing known V and J gene segments and tags...')

        # Look for tag and V/J fasta files: if cannot find locally, sources from GitHub repositories
        if tags == "original":
          if os.path.isfile("human_TRAV_region.fasta"):
            handle = open("human_TRAV_region.fasta", "rU")
          else:
            handle = urllib2.urlopen("https://raw.githubusercontent.com/uclinfectionimmunity/Decombinator/master/human_TRAV_region.fasta")
        elif tags == "extended":
          if os.path.isfile("exthuman_TRAV_region.fasta"):  
            handle = open("exthuman_TRAV_region.fasta", "rU")
          else:
            handle = urllib2.urlopen("https://raw.githubusercontent.com/JamieHeather/tcr-analysis/master/exthuman_TRAV_region.fasta")
        v_genes = list(SeqIO.parse(handle, "fasta"))
        handle.close()

        if tags == "original":
          if os.path.isfile("human_TRAJ_region.fasta"):
            handle = open("human_TRAJ_region.fasta", "rU")
          else:
            handle =    urllib2.urlopen("https://raw.githubusercontent.com/uclinfectionimmunity/Decombinator/master/human_TRAJ_region.fasta")    
        elif tags == "extended":
          if os.path.isfile("exthuman_TRAJ_region.fasta"):  
            handle = open("exthuman_TRAJ_region.fasta", "rU")
          else:
            handle = urllib2.urlopen("https://raw.githubusercontent.com/JamieHeather/tcr-analysis/master/exthuman_TRAJ_region.fasta")    
        j_genes = list(SeqIO.parse(handle, "fasta"))
        handle.close()

        v_regions = []
        for j in range(0, len(v_genes)):
            v_regions.append(string.upper(v_genes[j].seq))

        j_regions = []
        for j in range(0, len(j_genes)):
            j_regions.append(string.upper(j_genes[j].seq))

        ##############
        ## Build keyword tries of V and J tags for fast assignment
        if tags == "original":
          if os.path.isfile("tags_trav.txt") and os.path.isfile("tags_traj.txt"):
            v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(open("tags_trav.txt", "rU"), v_half_split)
            j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(open("tags_traj.txt", "rU"), j_half_split)
          else:
            v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(urllib2.urlopen("https://raw.githubusercontent.com/uclinfectionimmunity/Decombinator/master/humantags_trav.txt"), v_half_split)
            j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(urllib2.urlopen("https://raw.githubusercontent.com/uclinfectionimmunity/Decombinator/master/humantags_traj.txt"), j_half_split)               
        elif tags == "extended":
          if os.path.isfile("exttags_trav.txt") and os.path.isfile("exttags_traj.txt"):  
            v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(open("exttags_trav.txt", "rU"), v_half_split)
            j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(open("exttags_traj.txt", "rU"), j_half_split)
          else:
            v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(urllib2.urlopen("https://raw.githubusercontent.com/JamieHeather/tcr-analysis/master/exttags_trav.txt"), v_half_split)
            j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(urllib2.urlopen("https://raw.githubusercontent.com/JamieHeather/tcr-analysis/master/exttags_traj.txt"), j_half_split) 

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

        ###############
        ## Initialise variables
        assigned_count = 0 # this will just increase by one every time we correctly assign a seq read with all desired variables
        seq_count = 0 # this will simply track the number of sequences analysed in file
        t0 = time() # Begin timer

        ###############
        ## Open .txt file created at the start of analysis
        stemplate = Template('$v $j $del_v $del_j $nt_insert $seqid $tcr_seq $tcr_qual $barcode $barqual')      
 # Creates stemplate, a holder, for f. Each line will have the 5 variables separated by a space
 

        ###############
        ## Begin analysing sequences

        for i in range(len(Sequence_Reads)):
            
            print 'Importing sequences from', Sequence_Reads[i],' and assigning V and J regions...'
            
            if sys.stdin.isatty() == True:
              handle = open(Sequence_Reads[i], "rU") 
            else:
              handle = sys.stdin
            
            for record in SeqIO.parse(handle, "fastq"):
              
                barcode_seq = str(record.seq)[:12]         
                barcode_qual = str(record.format("fastq")).split("\n")[3][:12]          
                
                if 'N' in str(record.seq) and omitN==True:
                    Nseqs += 1
                else:
                    found_seq_match = 0
                    found_v_match = 0
                    found_j_match = 0
                    seq_count += 1

                    v_seq_start = 0     
                    j_seq_end = 0                          
                    
                    
                    hold_v = v_key.findall(str(record.seq))
                    hold_j = j_key.findall(str(record.seq))

                    if hold_v:
                        v_match = v_seqs.index(hold_v[0][0]) # Assigns V
                        temp_end_v = hold_v[0][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                        
                        v_seq_start = hold_v[0][1]      
                        
                        
                        if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ): # If the number of deletions has been found
                            [ end_v, deletions_v] = get_v_deletions( record.seq, v_match, temp_end_v, v_regions )
                    else:
                        found_v_match = 0
                        hold_v1 = half1_v_key.findall(str(record.seq))
                        hold_v2 = half2_v_key.findall(str(record.seq))
                        for i in range(len(hold_v1)):
                            indices = [y for y, x in enumerate(half1_v_seqs) if x == hold_v1[i][0] ]
                            for k in indices:
                                if len(v_seqs[k]) == len(str(record.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[half1_v_seqs.index(hold_v1[i][0])])]):
                                    if lev.hamming( v_seqs[k], str(record.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[k])] ) <= 1:
                                        v_match = k
                                        temp_end_v = hold_v1[i][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                                        found_v_match += 1
                                        
                                        v_seq_start = hold_v1[i][1]      

                                        
                                        
                        for i in range(len(hold_v2)):
                            indices = [y for y, x in enumerate(half2_v_seqs) if x == hold_v2[i][0] ]
                            for k in indices:
                                if len(v_seqs[k]) == len(str(record.seq)[hold_v2[i][1]-v_half_split:hold_v2[i][1]-v_half_split+len(v_seqs[half2_v_seqs.index(hold_v2[i][0])])]):
                                    if lev.hamming( v_seqs[k], str(record.seq)[hold_v2[i][1]-v_half_split:hold_v2[i][1]+len(v_seqs[k])-v_half_split] ) <= 1:
                                        v_match = k
                                        temp_end_v = hold_v2[i][1] + jump_to_end_v[v_match] - v_half_split - 1 # Finds where the end of a full V would be
                                        found_v_match += 1
                                        
                                        v_seq_start = hold_v2[i][1] - v_half_split      
                                        

                    if hold_j:
                        j_match = j_seqs.index(hold_j[0][0]) # Assigns J
                        temp_start_j = hold_j[0][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
                        
                        j_seq_end = hold_j[0][1] + len(hold_j[0][0])      
                        
                        if get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                            [ start_j, deletions_j] = get_j_deletions( record.seq, j_match, temp_start_j, j_regions )
                    else:
                        found_j_match = 0
                        hold_j1 = half1_j_key.findall(str(record.seq))
                        hold_j2 = half2_j_key.findall(str(record.seq))
                        for i in range(len(hold_j1)):
                            indices = [y for y, x in enumerate(half1_j_seqs) if x == hold_j1[i][0] ]
                            for k in indices:
                                if len(j_seqs[k]) == len(str(record.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[half1_j_seqs.index(hold_j1[i][0])])]):
                                    if lev.hamming( j_seqs[k], str(record.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[k])] ) <= 1:
                                        j_match = k
                                        temp_start_j = hold_j1[i][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
                                        found_j_match += 1

                                        j_seq_end = hold_j1[i][1] + len(hold_j1[i][0]) + j_half_split                                              
                                        
                        for i in range(len(hold_j2)):
                            indices = [y for y, x in enumerate(half2_j_seqs) if x == hold_j2[i][0] ]
                            for k in indices:
                                if len(j_seqs[k]) == len(str(record.seq)[hold_j2[i][1]-j_half_split:hold_j2[i][1]-j_half_split+len(j_seqs[half2_j_seqs.index(hold_j2[i][0])])]):
                                    if lev.hamming( j_seqs[k], str(record.seq)[hold_j2[i][1]-j_half_split:hold_j2[i][1]+len(j_seqs[k])-j_half_split] ) <= 1:
                                        j_match = k
                                        temp_start_j = hold_j2[i][1] - jump_to_start_j[j_match] - j_half_split # Finds where the start of a full J would be
                                        found_j_match += 1
                                        
                                        j_seq_end = hold_j2[i][1] + len(hold_j2[i][0])                                                

                    TCRseq = str(record.seq[v_seq_start:j_seq_end])       
                    TCRqual = str(record.format("fastq")).split("\n")[3][v_seq_start:j_seq_end]         

                    if hold_v and hold_j:
                        if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                            f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )  
                            
                            if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                imposs_deletion += 1                    
                            elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                tag_overlap += 1                     
                            else:        
                                print >> results_file, f_seq          
                                assigned_count += 1                   
                                found_seq_match = 1                         
                                    
                    elif hold_v and found_j_match == 1:
                        if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                            [ start_j, deletions_j] = get_j_deletions( record.seq, j_match, temp_start_j, j_regions )
                            f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )  
                            
                            if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                imposs_deletion += 1                    
                            elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                tag_overlap += 1                     
                            else:        
                                print >> results_file, f_seq          
                                assigned_count += 1                   
                                found_seq_match = 1                         
                            
                    elif found_v_match == 1 and hold_j:
                        if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                            [ end_v, deletions_v] = get_v_deletions( record.seq, v_match, temp_end_v, v_regions )
                            f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )  
                                                        
                            if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                imposs_deletion += 1                    
                            elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                tag_overlap += 1                     
                            else:        
                                print >> results_file, f_seq          
                                assigned_count += 1                   
                                found_seq_match = 1                         

                    elif found_v_match == 1 and found_j_match == 1:
                        if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                            [ end_v, deletions_v] = get_v_deletions( record.seq, v_match, temp_end_v, v_regions )
                            [ start_j, deletions_j] = get_j_deletions( record.seq, j_match, temp_start_j, j_regions )
                            f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )  
                            
                            if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                imposs_deletion += 1                    
                            elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                tag_overlap += 1                     
                            else:        
                                print >> results_file, f_seq          
                                assigned_count += 1                   
                                found_seq_match = 1                         
                            
                    if found_seq_match == 0 and with_reverse_complement_search == True:
                        
                        #####################
                        # REVERSE COMPLEMENT
                        #####################

                        record_reverse = record.reverse_complement()
                        hold_v = v_key.findall(str(record_reverse.seq))
                        hold_j = j_key.findall(str(record_reverse.seq))

                        if hold_v:
                            v_match = v_seqs.index(hold_v[0][0]) # Assigns V
                            temp_end_v = hold_v[0][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                            
                            v_seq_start = hold_v[0][1]      
                            
                            if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ): # If the number of deletions has been found
                                [ end_v, deletions_v] = get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions )
                        else:
                            found_v_match = 0
                            hold_v1 = half1_v_key.findall(str(record_reverse.seq))
                            hold_v2 = half2_v_key.findall(str(record_reverse.seq))
                            for i in range(len(hold_v1)):
                                indices = [y for y, x in enumerate(half1_v_seqs) if x == hold_v1[i][0] ]
                                for k in indices:
                                    if len(v_seqs[k]) == len(str(record_reverse.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[half1_v_seqs.index(hold_v1[i][0])])]):
                                        if lev.hamming( v_seqs[k], str(record_reverse.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[k])] ) <= 1:
                                            v_match = k
                                            temp_end_v = hold_v1[i][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                                            found_v_match += 1
                                            
                                            v_seq_start = hold_v1[i][1]      
                                            
                            for i in range(len(hold_v2)):
                                indices = [y for y, x in enumerate(half2_v_seqs) if x == hold_v2[i][0] ]
                                for k in indices:
                                    if len(v_seqs[k]) == len(str(record_reverse.seq)[hold_v2[i][1]-v_half_split:hold_v2[i][1]-v_half_split+len(v_seqs[half2_v_seqs.index(hold_v2[i][0])])]):
                                        if lev.hamming( v_seqs[k], str(record_reverse.seq)[hold_v2[i][1]-v_half_split:hold_v2[i][1]+len(v_seqs[k])-v_half_split] ) <= 1:
                                            v_match = k
                                            temp_end_v = hold_v2[i][1] + jump_to_end_v[v_match] - v_half_split - 1 # Finds where the end of a full V would be
                                            found_v_match += 1
                                            
                                            v_seq_start = hold_v2[i][1] - v_half_split                                                  

                        if hold_j:
                            j_match = j_seqs.index(hold_j[0][0]) # Assigns J
                            temp_start_j = hold_j[0][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
                            
                            j_seq_end = hold_j[0][1] + len(hold_j[0][0])      
                            
                            if get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                                [ start_j, deletions_j] = get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions )
                        else:
                            found_j_match = 0
                            hold_j1 = half1_j_key.findall(str(record_reverse.seq))
                            hold_j2 = half2_j_key.findall(str(record_reverse.seq))
                            for i in range(len(hold_j1)):
                                indices = [y for y, x in enumerate(half1_j_seqs) if x == hold_j1[i][0] ]
                                for k in indices:
                                    if len(j_seqs[k]) == len(str(record_reverse.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[half1_j_seqs.index(hold_j1[i][0])])]):
                                        if lev.hamming( j_seqs[k], str(record_reverse.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[k])] ) <= 1:
                                            j_match = k
                                            temp_start_j = hold_j1[i][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
                                            found_j_match += 1
                                            
                                            j_seq_end = hold_j1[i][1] + len(hold_j1[i][0]) + j_half_split      
                                            
                            for i in range(len(hold_j2)):
                                indices = [y for y, x in enumerate(half2_j_seqs) if x == hold_j2[i][0] ]
                                for k in indices:
                                    if len(j_seqs[k]) == len(str(record_reverse.seq)[hold_j2[i][1]-j_half_split:hold_j2[i][1]-j_half_split+len(j_seqs[half2_j_seqs.index(hold_j2[i][0])])]):
                                        if lev.hamming( j_seqs[k], str(record_reverse.seq)[hold_j2[i][1]-j_half_split:hold_j2[i][1]+len(j_seqs[k])-j_half_split] ) <= 1:
                                            j_match = k
                                            temp_start_j = hold_j2[i][1] - jump_to_start_j[j_match] - j_half_split # Finds where the start of a full J would be
                                            found_j_match += 1

                                            j_seq_end = hold_j2[i][1] + len(hold_j2[i][0])        

                        # Records the inter-tag sequence and quality for downstream error-correction
                        TCRseq = str(record_reverse.seq[v_seq_start:j_seq_end])       
                        TCRqual = str(record_reverse.format("fastq")).split("\n")[3][v_seq_start:j_seq_end]     
                        

                        if hold_v and hold_j:
                            if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                                f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record_reverse.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )      
                                
                                if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                    imposs_deletion += 1                    
                                elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                    tag_overlap += 1                     
                                else:        
                                    print >> results_file, f_seq          
                                    assigned_count += 1                   
                                    found_seq_match = 1                                                        

                        elif hold_v and found_j_match == 1:
                            if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                                [ start_j, deletions_j] = get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions )
                                f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record_reverse.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )      

                                if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                    imposs_deletion += 1                    
                                elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                    tag_overlap += 1                     
                                else:        
                                    print >> results_file, f_seq          
                                    assigned_count += 1                   
                                    found_seq_match = 1                                                        

                        elif found_v_match == 1 and hold_j:
                            if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                                [ end_v, deletions_v] = get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions )
                                f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record_reverse.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )      

                                if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                    imposs_deletion += 1                    
                                elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                    tag_overlap += 1                     
                                else:        
                                    print >> results_file, f_seq          
                                    assigned_count += 1                   
                                    found_seq_match = 1                        
                                    
                        elif found_v_match == 1 and found_j_match == 1:
                            if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                                [ end_v, deletions_v] = get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions )
                                [ start_j, deletions_j] = get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions )
                                f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record_reverse.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )      

                                if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                    imposs_deletion += 1                    
                                elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                    tag_overlap += 1                     
                                else:        
                                    print >> results_file, f_seq          
                                    assigned_count += 1                   
                                    found_seq_match = 1                        
                                    
            handle.close()
        results_file.close()
    
    elif chain=="beta":
        
         # Do not change - V tags are split at position 10, J at position 6 (for original tags)
          # and at position 10 for extended tags.
        if tags == "original":
          v_half_split, j_half_split = [10,6]
        elif tags == "extended":
          v_half_split, j_half_split = [10,10] 
        ################

        print 'Commencing analysis on a total of', len(Sequence_Reads), 'file(s)'

        ################
        print ('Importing known V, D and J gene segments and tags...')

        # Look for tag and V/J fasta files: if cannot find locally, sources from GitHub repositories
        if tags == "original":
          if os.path.isfile("human_TRBV_region.fasta"):
            handle = open("human_TRBV_region.fasta", "rU")
          else:
            handle = urllib2.urlopen("https://raw.githubusercontent.com/uclinfectionimmunity/Decombinator/master/human_TRBV_region.fasta")
        elif tags == "extended":
          if os.path.isfile("exthuman_TRBV_region.fasta"):        
            handle = open("exthuman_TRBV_region.fasta", "rU")
          else:
            handle = urllib2.urlopen("https://raw.githubusercontent.com/JamieHeather/tcr-analysis/master/exthuman_TRBV_region.fasta")
        v_genes = list(SeqIO.parse(handle, "fasta"))
        handle.close()

        if tags == "original":
          if os.path.isfile("human_TRBJ_region.fasta"):   
            handle = open("human_TRBJ_region.fasta", "rU")
          else:
            handle = urllib2.urlopen("https://raw.githubusercontent.com/uclinfectionimmunity/Decombinator/master/human_TRBJ_region.fasta")    
        elif tags == "extended":        
          if os.path.isfile("exthuman_TRBJ_region.fasta"):  
            handle = open("exthuman_TRBJ_region.fasta", "rU")
          else:
            handle = urllib2.urlopen("https://raw.githubusercontent.com/JamieHeather/tcr-analysis/master/exthuman_TRBJ_region.fasta")    
        j_genes = list(SeqIO.parse(handle, "fasta"))
        handle.close()

        v_regions = []
        for j in range(0, len(v_genes)):
            v_regions.append(string.upper(v_genes[j].seq))

        j_regions = []
        for j in range(0, len(j_genes)):
            j_regions.append(string.upper(j_genes[j].seq))

        ##############
        ## Build keyword tries of V and J tags for fast assignment
        if tags == "original":
          if os.path.isfile("tags_trbv.txt") and os.path.isfile("tags_trbj.txt"):
            v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(open("tags_trbv.txt", "rU"), v_half_split)
            j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(open("tags_trbj.txt", "rU"), j_half_split)
          else:
            v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(urllib2.urlopen("https://raw.githubusercontent.com/uclinfectionimmunity/Decombinator/master/humantags_trbv.txt"), v_half_split)
            j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(urllib2.urlopen("https://raw.githubusercontent.com/uclinfectionimmunity/Decombinator/master/humantags_trbj.txt"), j_half_split)        
        elif tags == "extended":
          if os.path.isfile("exttags_trbv.txt") and os.path.isfile("exttags_trbj.txt"):  
            v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(open("exttags_trbv.txt", "rU"), v_half_split)
            j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(open("exttags_trbj.txt", "rU"), j_half_split)
          else:
            v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(urllib2.urlopen("https://raw.githubusercontent.com/JamieHeather/tcr-analysis/master/exttags_trbv.txt"), v_half_split)
            j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(urllib2.urlopen("https://raw.githubusercontent.com/JamieHeather/tcr-analysis/master/exttags_trbj.txt"), j_half_split)       
  
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

        ###############
        ## Initialise variables
        assigned_count = 0 # this will just increase by one every time we correctly assign a seq read with all desired variables
        seq_count = 0 # this will simply track the number of sequences analysed in file
        Nseqs = 0 # Number of raw reads containing N nucleotides
        t0 = time() # Begin timer

        ###############
        ## Open .txt file created at the start of analysis
        stemplate = Template('$v $j $del_v $del_j $nt_insert $seqid $tcr_seq $tcr_qual $barcode $barqual')      
 # Creates stemplate, a holder, for f. Each line will have the 5 variables separated by a space

        ###############
        ## Begin analysing sequences

        for i in range(len(Sequence_Reads)):
            
            print 'Importing sequences from', Sequence_Reads[i],' and assigning V and J regions...'
            
            if sys.stdin.isatty() == True:
              handle = open(Sequence_Reads[i], "rU") 
            else:
              handle = sys.stdin
            
            for record in SeqIO.parse(handle, "fastq"):
              
                barcode_seq = str(record.seq)[0:12]         
                barcode_qual = str(record.format("fastq")).split("\n")[3][:12]          
                
                if 'N' in str(record.seq) and omitN==True:
                    Nseqs += 1
                else:
                    found_seq_match = 0
                    found_v_match = 0
                    found_j_match = 0
                    seq_count += 1
                                      
                    v_seq_start = 0     
                    j_seq_end = 0       
                                       
                    hold_v = v_key.findall(str(record.seq))
                    hold_j = j_key.findall(str(record.seq))

                    if hold_v:
                        v_match = v_seqs.index(hold_v[0][0]) # Assigns V
                        temp_end_v = hold_v[0][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                        
                        v_seq_start = hold_v[0][1]      
                        
                        if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ): # If the number of deletions has been found
                            [ end_v, deletions_v] = get_v_deletions( record.seq, v_match, temp_end_v, v_regions )
                    else:
                        found_v_match = 0
                        hold_v1 = half1_v_key.findall(str(record.seq))
                        hold_v2 = half2_v_key.findall(str(record.seq))
                        for i in range(len(hold_v1)):
                            indices = [y for y, x in enumerate(half1_v_seqs) if x == hold_v1[i][0] ]
                            for k in indices:
                                if len(v_seqs[k]) == len(str(record.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[half1_v_seqs.index(hold_v1[i][0])])]):
                                    if lev.hamming( v_seqs[k], str(record.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[k])] ) <= 1:
                                        v_match = k
                                        temp_end_v = hold_v1[i][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                                        found_v_match += 1
                                        
                                        v_seq_start = hold_v1[i][1]                                              
                                        
                        for i in range(len(hold_v2)):
                            indices = [y for y, x in enumerate(half2_v_seqs) if x == hold_v2[i][0] ]
                            for k in indices:
                                if len(v_seqs[k]) == len(str(record.seq)[hold_v2[i][1]-v_half_split:hold_v2[i][1]+len(v_seqs[half2_v_seqs.index(hold_v2[i][0])])-v_half_split]):
                                    if lev.hamming( v_seqs[k], str(record.seq)[hold_v2[i][1]-v_half_split:hold_v2[i][1]+len(v_seqs[k])-v_half_split] ) <= 1:
                                        v_match = k
                                        temp_end_v = hold_v2[i][1] + jump_to_end_v[v_match] - v_half_split - 1 # Finds where the end of a full V would be
                                        found_v_match += 1

                                        v_seq_start = hold_v2[i][1] - v_half_split      
                                        

                    if hold_j:
                        j_match = j_seqs.index(hold_j[0][0]) # Assigns J
                        temp_start_j = hold_j[0][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
                        
                        j_seq_end = hold_j[0][1] + len(hold_j[0][0])      
                        
                        if get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                            [ start_j, deletions_j] = get_j_deletions( record.seq, j_match, temp_start_j, j_regions )
                    else:
                        found_j_match = 0
                        hold_j1 = half1_j_key.findall(str(record.seq))
                        hold_j2 = half2_j_key.findall(str(record.seq))
                        for i in range(len(hold_j1)):
                            indices = [y for y, x in enumerate(half1_j_seqs) if x == hold_j1[i][0] ]
                            for k in indices:
                                if len(j_seqs[k]) == len(str(record.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[half1_j_seqs.index(hold_j1[i][0])])]):
                                    if lev.hamming( j_seqs[k], str(record.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[k])] ) <= 1:
                                        j_match = k
                                        temp_start_j = hold_j1[i][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
                                        found_j_match += 1
                                        
                                        ## need to account for the fact that there are different length TRBJ tags! 
                                        #if j_match == 7:                                         
                                            #j_seq_end = hold_j1[i][1] + len(hold_j1[i][0]) + j_half_split + 2    
                                        #elif j_match == 12:                    
                                            #j_seq_end = hold_j1[i][1] + len(hold_j1[i][0]) + j_half_split + 3     
                                        #else:                    
                                        j_seq_end = hold_j1[i][1] + len(hold_j1[i][0]) + j_half_split      
                                        
                        for i in range(len(hold_j2)):
                            indices = [y for y, x in enumerate(half2_j_seqs) if x == hold_j2[i][0] ]
                            for k in indices:
                                if len(j_seqs[k]) == len(str(record.seq)[hold_j2[i][1]-j_half_split:hold_j2[i][1]+len(j_seqs[half2_j_seqs.index(hold_j2[i][0])])-j_half_split]):
                                    if lev.hamming( j_seqs[k], str(record.seq)[hold_j2[i][1]-j_half_split:hold_j2[i][1]+len(j_seqs[k])-j_half_split] ) <= 1:
                                        j_match = k
                                        temp_start_j = hold_j2[i][1] - jump_to_start_j[j_match] - j_half_split # Finds where the start of a full J would be
                                        found_j_match += 1
                                        
                                        j_seq_end = hold_j2[i][1] + len(hold_j2[i][0])        
                                        

                    TCRseq = str(record.seq[v_seq_start:j_seq_end])       
                    TCRqual = str(record.format("fastq")).split("\n")[3][v_seq_start:j_seq_end]                                 
                                                            

                    if hold_v and hold_j:
                        if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                            f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )                              

                            if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                imposs_deletion += 1                    
                            elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                tag_overlap += 1                     
                            else:        
                                print >> results_file, f_seq          
                                assigned_count += 1                   
                                found_seq_match = 1                         

                    elif hold_v and found_j_match == 1:
                        if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                            [ start_j, deletions_j] = get_j_deletions( record.seq, j_match, temp_start_j, j_regions )
                            f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )  
                            
                            if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                imposs_deletion += 1                    
                            elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                tag_overlap += 1                     
                            else:        
                                print >> results_file, f_seq          
                                assigned_count += 1                   
                                found_seq_match = 1                         

                    elif found_v_match == 1 and hold_j:
                        if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                            [ end_v, deletions_v] = get_v_deletions( record.seq, v_match, temp_end_v, v_regions )
                            f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )  
                            
                            if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                imposs_deletion += 1                    
                            elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                tag_overlap += 1                     
                            else:        
                                print >> results_file, f_seq          
                                assigned_count += 1                   
                                found_seq_match = 1                         

                    elif found_v_match == 1 and found_j_match == 1:
                        if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                            [ end_v, deletions_v] = get_v_deletions( record.seq, v_match, temp_end_v, v_regions )
                            [ start_j, deletions_j] = get_j_deletions( record.seq, j_match, temp_start_j, j_regions )
                            f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )  
                            
                            if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                imposs_deletion += 1                    
                            elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                tag_overlap += 1                     
                            else:        
                                print >> results_file, f_seq          
                                assigned_count += 1                   
                                found_seq_match = 1                         


                    if found_seq_match == 0 and with_reverse_complement_search == True:
                        
                        #####################
                        # REVERSE COMPLEMENT
                        #####################

                        record_reverse = record.reverse_complement()
                        hold_v = v_key.findall(str(record_reverse.seq))
                        hold_j = j_key.findall(str(record_reverse.seq))

                        if hold_v:
                            v_match = v_seqs.index(hold_v[0][0]) # Assigns V
                            temp_end_v = hold_v[0][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                            
                            v_seq_start = hold_v[0][1]      
                            
                            if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ): # If the number of deletions has been found
                                [ end_v, deletions_v] = get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions )
                        else:
                            found_v_match = 0
                            hold_v1 = half1_v_key.findall(str(record_reverse.seq))
                            hold_v2 = half2_v_key.findall(str(record_reverse.seq))
                            for i in range(len(hold_v1)):
                                indices = [y for y, x in enumerate(half1_v_seqs) if x == hold_v1[i][0] ]
                                for k in indices:
                                    if len(v_seqs[k]) == len(str(record_reverse.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[half1_v_seqs.index(hold_v1[i][0])])]):
                                        if lev.hamming( v_seqs[k], str(record_reverse.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[k])] ) <= 1:
                                            v_match = k
                                            temp_end_v = hold_v1[i][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                                            found_v_match += 1
                                            
                                            v_seq_start = hold_v1[i][1]      
                                            
                            for i in range(len(hold_v2)):
                                indices = [y for y, x in enumerate(half2_v_seqs) if x == hold_v2[i][0] ]
                                for k in indices:
                                    if len(v_seqs[k]) == len(str(record_reverse.seq)[hold_v2[i][1]-v_half_split:hold_v2[i][1]+len(v_seqs[half2_v_seqs.index(hold_v2[i][0])])-v_half_split]):
                                        if lev.hamming( v_seqs[k], str(record_reverse.seq)[hold_v2[i][1]-v_half_split:hold_v2[i][1]+len(v_seqs[k])-v_half_split] ) <= 1:
                                            v_match = k
                                            temp_end_v = hold_v2[i][1] + jump_to_end_v[v_match] - v_half_split - 1 # Finds where the end of a full V would be
                                            found_v_match += 1
                                            
                                            v_seq_start = hold_v2[i][1] - v_half_split                                                  

                        if hold_j:
                            j_match = j_seqs.index(hold_j[0][0]) # Assigns J
                            temp_start_j = hold_j[0][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
                            
                            j_seq_end = hold_j[0][1] + len(hold_j[0][0])      
                            
                            if get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                                [ start_j, deletions_j] = get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions )
                        else:
                            found_j_match = 0
                            hold_j1 = half1_j_key.findall(str(record_reverse.seq))
                            hold_j2 = half2_j_key.findall(str(record_reverse.seq))
                            for i in range(len(hold_j1)):
                                indices = [y for y, x in enumerate(half1_j_seqs) if x == hold_j1[i][0] ]
                                for k in indices:
                                    if len(j_seqs[k]) == len(str(record_reverse.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[half1_j_seqs.index(hold_j1[i][0])])]):
                                        if lev.hamming( j_seqs[k], str(record_reverse.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[k])] ) <= 1:
                                            j_match = k
                                            temp_start_j = hold_j1[i][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
                                            found_j_match += 1
                                            
                                            # need to account for the fact that there are different length TRBJ tags! 
                                            #if j_match == 7:                                         
                                                #j_seq_end = hold_j1[i][1] + len(hold_j1[i][0]) + j_half_split + 2    
                                            #elif j_match == 12:                    
                                                #j_seq_end = hold_j1[i][1] + len(hold_j1[i][0]) + j_half_split + 3     
                                            #else:                    
                                            j_seq_end = hold_j1[i][1] + len(hold_j1[i][0]) + j_half_split      
                                                                                            
                            for i in range(len(hold_j2)):
                                indices = [y for y, x in enumerate(half2_j_seqs) if x == hold_j2[i][0] ]
                                for k in indices:
                                    if len(j_seqs[k]) == len(str(record_reverse.seq)[hold_j2[i][1]-j_half_split:hold_j2[i][1]+len(j_seqs[half2_j_seqs.index(hold_j2[i][0])])-j_half_split]):
                                        if lev.hamming( j_seqs[k], str(record_reverse.seq)[hold_j2[i][1]-j_half_split:hold_j2[i][1]+len(j_seqs[k])-j_half_split] ) <= 1:
                                            j_match = k
                                            temp_start_j = hold_j2[i][1] - jump_to_start_j[j_match] - j_half_split # Finds where the start of a full J would be
                                            found_j_match += 1

                                            j_seq_end = hold_j2[i][1] + len(hold_j2[i][0])        

                        # Records the inter-tag sequence and quality for downstream error-correction                                            
                        TCRseq = str(record_reverse.seq[v_seq_start:j_seq_end])       
                        TCRqual = str(record_reverse.format("fastq")).split("\n")[3][v_seq_start:j_seq_end]     


                        if hold_v and hold_j:
                            if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                                f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record_reverse.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )      
                                

                                if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                    imposs_deletion += 1        
                                elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                    tag_overlap += 1        
                                else:        
                                    print >> results_file, f_seq        
                                    assigned_count += 1        
                                    found_seq_match = 1               
                                    
                        elif hold_v and found_j_match == 1:
                            if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                                [ start_j, deletions_j] = get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions )
                                f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record_reverse.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )      
     
                                if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                    imposs_deletion += 1        
                                elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                    tag_overlap += 1        
                                else:        
                                    print >> results_file, f_seq        
                                    assigned_count += 1        
                                    found_seq_match = 1                   
     
                        elif found_v_match == 1 and hold_j:
                            if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                                [ end_v, deletions_v] = get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions )
                                f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record_reverse.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )      

                                if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                    imposs_deletion += 1                    
                                elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                    tag_overlap += 1                     
                                else:        
                                    print >> results_file, f_seq          
                                    assigned_count += 1                   
                                    found_seq_match = 1                                              
                                
                        elif found_v_match == 1 and found_j_match == 1:
                            if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                                [ end_v, deletions_v] = get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions )
                                [ start_j, deletions_j] = get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions )
                                f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record_reverse.seq[end_v+1:start_j])+str(','), seqid = str(record.id)+str(','), tcr_seq = TCRseq+str(','), tcr_qual = TCRqual+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )      

                                if deletions_v > (jump_to_end_v[v_match] - len(v_seqs[v_match])) or deletions_j > jump_to_start_j[j_match]: # Impossible deletion filter        
                                    imposs_deletion += 1        
                                elif ((temp_end_v - jump_to_end_v[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_j[j_match]): # overlapping tag filter        
                                    tag_overlap += 1        
                                else:        
                                    print >> results_file, f_seq        
                                    assigned_count += 1        
                                    found_seq_match = 1                                              
                                
            handle.close()
        results_file.close()
   
    timed = time() - t0
    print seq_count, 'sequences were analysed' 
    print assigned_count, 'sequences were successfully assigned'
    if omitN==True:
        print Nseqs, 'sequences contained ambiguous N nucleotides'
    print tag_overlap, "reads with overlapping tags discarded"                       
    print imposs_deletion, "reads with impossible numbers of deletions discarded"       
    print "Output to " + name_results + ".n12"
         
    print 'Time taken =', timed, 'seconds'

###################
# OTHER FUNCTIONS #
###################

def get_v_deletions( rc, v_match, temp_end_v, v_regions_cut ):
    # This function determines the number of V deletions in sequence rc
    # by comparing it to v_match, beginning by making comparisons at the
    # end of v_match and at position temp_end_v in rc.
    function_temp_end_v = temp_end_v
    pos = -1
    is_v_match = 0
    while is_v_match == 0 and 0 <= function_temp_end_v < len(rc):
        if str(v_regions_cut[v_match])[pos] == str(rc)[function_temp_end_v] and str(v_regions_cut[v_match])[pos-1] == str(rc)[function_temp_end_v-1] and str(v_regions_cut[v_match])[pos-2] == str(rc)[function_temp_end_v-2]:
            is_v_match = 1
            deletions_v = -pos - 1
            end_v = function_temp_end_v
        else:
            pos -= 1
            function_temp_end_v -= 1

    if is_v_match == 1:
        return [end_v, deletions_v]
    else:
        return []

def get_j_deletions( rc, j_match, temp_start_j, j_regions_cut ):
    # This function determines the number of J deletions in sequence rc
    # by comparing it to j_match, beginning by making comparisons at the
    # end of j_match and at position temp_end_j in rc.
    function_temp_start_j = temp_start_j
    pos = 0
    is_j_match = 0
    while is_j_match == 0 and 0 <= function_temp_start_j+2 < len(str(rc)):
        if str(j_regions_cut[j_match])[pos] == str(rc)[function_temp_start_j] and str(j_regions_cut[j_match])[pos+1] == str(rc)[function_temp_start_j+1] and str(j_regions_cut[j_match])[pos+2] == str(rc)[function_temp_start_j+2]:
            is_j_match = 1
            deletions_j = pos
            start_j = function_temp_start_j
        else:
            pos += 1
            function_temp_start_j += 1
            
    if is_j_match == 1:
        return [start_j, deletions_j]
    else:
        return []

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


##################
### TAKE INPUT ###
##################

if (len(sys.argv) <> 3):
  print "Please enter both a file-name and which chain to decombine (a/b)"
  print "E.g. python vDCR.py FILE.fq a"
  sys.exit()
else:
  filename = str(sys.argv[1])
  chaininput = str(sys.argv[2])
  if chaininput == "a":
    chain = "alpha"
  elif chaininput == "b":
    chain = "beta"
  else:
    print "Chain must be \'a\' for alpha or \'b\' for beta"
    sys.exit()

    
print 'Running vDCR on ' + filename + '...'

filestoanalyse = []


try:
    testopen = open(str(filename),"rU")
    testopen.close()
    correctpath = True
    filestoanalyse.append(filename)
except:
    print 'Cannot find the file you specified. Please try again'
    sys.exit()
        
        
if tags == "original":
  name_results = "vDCR_" + str(chain) + "_" + str(filename.split(".")[0]) 
elif tags == "extended":
  name_results = "vDCRe_" + str(chain) + "_" + str(filename.split(".")[0]) 
else:
  print "Name of tag set not recognised. Please edit code so that tags = either \'original\' or \'extended\'"
  sys.exit()
  
currentpath = os.getcwd()

try:
    analysis( filestoanalyse, name_results, str(chain), with_statistics=True, with_reverse_complement_search=True, omitN=False)
except Exception, e:
    print 'verbose Decombinator encountered an unexpected error while analysing your file.'
    print 'CRYPTIC PYTHON ERROR MESSAGE SAYS:'
    print traceback.format_exc()

seqs_found = 0
seqcheck = open(name_results+'.n12',"rU")  
for line in seqcheck:
    seqs_found+=1
seqcheck.close()

if seqs_found == 0:
    print 'Could not find any TCR',chain,'sequences in the specified file.'

sys.exit()
