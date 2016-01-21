# testn12.2.py
# 21st Jan 2016 JH UCL
# want to get an idea of many of our decombined reads have the the barcodes in the right place, based on their barcode field

from __future__ import division
import gzip
import sys
import collections as coll
import Levenshtein as lev
import regex

infile = sys.argv[1]


def getbarcode(bcseq):
    """
    getbarcode(bcseq):
    Given a barcode-region sequence, outputs the sequence of the do-docamer barcode.
    This barcode (theoretically) consists of the concatentation of the two random hexamer sequences contained in the ligation oligo.
    However errors in sequences and ligation oligo production can mean that the random nucleotides are not always at the expected position.
    This function uses the known sequence of the spacers (which bound each of the two N6s to their 5') to deduce the random sequences
    """
    spcr = "GTCGTGAT"
    
    if "N" in bcseq:    # ambiguous base-call check
      counts['N'] += 1
      return
    
    # define expected positions
    first = bc[0:8]
    second = bc[14:22]   
    
    # if both spacers present and correct, return sequence from expected sites
    if first == spcr and second == spcr:
      n1 = bcseq[8:14] 
      n2 = bcseq[22:28]
      counts['getbarcode_pass_exactmatch'] += 1
      return n1+n2
    
    # otherwise look throughout the entire sequence for the presence of two spacers
    else:
      
      # pad the 5' of the sequence to allow for frame-shifts that result in incomplete primary spacer
      bcseq = "XXXX" + bcseq 
      
      # search for all instances of the spacer, allowing for 2 substitutions OR (1 deletion or 1 insertion)
      err_spcrs = regex.findall("(GTCGTGAT){2i+2d+1s<=2}", bcseq) 
      positions = [bcseq.find(x) for x in err_spcrs]
      lens = [len(x) for x in err_spcrs]
      
      if len(positions) <> 2:
        counts['getbarcode_fail_not2spacersfound'] += 1
        return
      else:
        # if only two matches, first random seq runs from end of first spacer to start of next
          # second random seq just
        n1 = bcseq[ positions[0]+lens[0] : positions[1] ]
        n2 = bcseq[ positions[1]+lens[1] : positions[1]+lens[1]+6 ]
        
        if len(n1) < 4:
          counts['getbarcode_fail_n1tooshort'] += 1
          return
        elif len(n1) > 8:
          counts['getbarcode_fail_n1toolong'] += 1
          
        concat = n1 + n2
        
        # need to output a 12-mer for downstream code
        if len(concat) == 12:
          counts['getbarcode_pass_fuzzymatch_rightlen'] += 1
          return concat
        elif len(concat) < 12:
          # if too short, pad sequence with 3' "S" characters (for "short")
          counts['getbarcode_pass_fuzzymatch_short'] += 1
          difflen = 12 - len(concat)
          return concat + ("S" * difflen)
        elif len(concat) > 12: 
          # if too long, chop back to 11 and append with a 3' "L" (for "long") 
            # this prevents clashing with a pre-existing shorter barcode
          counts['getbarcode_pass_fuzzymatch_long'] += 1
          return concat[:11] + "L"
        else:
          counts['getbarcode_fail_unknown'] += 1
          return
      
      
if infile.endswith('.gz'):
  opener = gzip.open
else:
  opener = open    

counts = coll.Counter()
nss = coll.Counter()

with opener(infile) as inf:
  
  for line in inf:
    
    counts['count'] += 1
     
    splt = line.rstrip().split(", ")
    bc = splt[8]
    
    ns = getbarcode(bc)
    
    if ns:
      #print ns
      nss[ns] += 1
      counts['success'] += 1
    else:
      counts['fail'] += 1
      
print nss.most_common(10)
print counts
