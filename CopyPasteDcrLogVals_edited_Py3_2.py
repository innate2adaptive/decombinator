# Last updated : July 2017

# MI



##### Pseudocode #####



# read in Decombinator log files

# for each file

 # store Outputfile - this will be what we look for in the Coll log file

 # store the values we want from the Decombinator log file

# read in Collapsing log files

 # find corresponding Outputfile name in the log file

 # store the values we want from the Collapsing log file



# output:

 # sample,NumberReadsInput,NumberReadsDecombined,PercentReadsDecombined,UniqueDCRsPassingFilters,TotalDCRsPassingFilters,PercentDCRPassingFilters(withbarcode),UniqueDCRsPostCollapsing,TotalDCRsPostCollapsing,PercentUniqueDCRsKept,PercentTotalDCRsKept,AverageInputTCRAbundance,AverageOutputTCRAbundance,AverageRNAduplication

 

# output comes from:

 # field 1 = Decombinator log, field 2-4 = Decombinator log, field 5-14 = Collapsing log



##### How to run #####



# python CopyPasteDcrLogVals.py full/path/to/Logs/folder/ outfile.csv



# NB 

# assumes that all Log files are in one Log folder

# prints the output file to whatever directory you are in

# may require slight modification on line 94, depending on character separator 

#  for grabbing filename in original log file

 

##### Py packages #####



from os import listdir

from os.path import isfile, join
import sys
 
##### Start #####

pathToLogs = str(sys.argv[1]) # full path to Logs folder

onlyfiles = [f for f in listdir(pathToLogs) if isfile(join(pathToLogs, f))]

sampleNam = []

# get sample names from Decombinator log file
# then get values 

for f in onlyfiles:

    if 'Decombinator' in f:

        with open(pathToLogs+f, 'r') as f:
            lines = f.read().splitlines()

            for l in lines:

                if 'OutputFile' in l:
                    spl = l.split(',')
                    sampleNam.append(spl[1])                  

fields = ["sample",
          "NumberReadsInput",
          "NumberReadsDecombined",
          "PercentReadsDecombined",
          "UniqueDCRsPassingFilters",
          "TotalDCRsPassingFilters",
          "PercentDCRPassingFilters(withbarcode)",
          "UniqueDCRsPostCollapsing",
          "TotalDCRsPostCollapsing",
          "PercentUniqueDCRsKept",
          "PercentTotalDCRsKept",
          "AverageInputTCRAbundance",
          "AverageOutputTCRAbundance",
          "AverageRNAduplication"]

out = []

for i in sampleNam:

    string = [i]

    for j in onlyfiles:
        if 'Decombinator' in j:
            with open(pathToLogs+j, 'r') as inf:
                lines = inf.read().splitlines()

                for idx, l in enumerate(lines):

                    if 'OutputFile' in l:
                        nam = l.split(',')[1]

                        if nam == i:

                            for item in lines[idx:]:
                                for patt in fields[1:4]:
                                    if patt in item:
                                        val = item.split(',')[1]
                                        string.append(val)

    for j in onlyfiles:

        if 'Collapsing' in j:

            with open(pathToLogs+j, 'r') as inf:
                lines = inf.read().splitlines()

                for idx, l in enumerate(lines):
                    if 'InputFile' in l:
                        nam = l.split('/')[-1] # modify here for either
                                              # comma or slash

                        if nam == i:
                            for item in lines[idx:]:
                                for patt in fields[4:]:
                                    if patt in item:
                                        val = item.split(',')[1]
                                        string.append(val)

    # modify filename fun

    for char in ['.', '-']:

        if char in string[0]:
            string[0] = string[0].replace(char, '_')

    spl = string[0].split('_')[2:-2]
    a_idx = [i for i, x, in enumerate(spl) if x == 'a']
    b_idx = [i for i, x, in enumerate(spl) if x == 'b']

    # print(spl, a_idx, b_idx)

    if string[0].startswith('dcr_alpha'):

        for i, j in zip(a_idx, b_idx):

            if i < j:

                new_nam = 'alpha_' + '_'.join(spl[:i+1])
                string[0] = new_nam

            if i > j:

                new_nam = 'alpha_' + '_'.join(spl[j+1:])
                string[0] = new_nam

        if not b_idx:

            new_nam = 'alpha_' + '_'.join(spl)
            string[0] = new_nam

    if string[0].startswith('dcr_beta'):
        
        for i, j in zip(a_idx, b_idx):

            if i < j:
                new_nam = 'beta_' + '_'.join(spl[i+1:])
                string[0] = new_nam

            if i > j:
                new_nam = 'beta_' + '_'.join(spl[:j+1])
                string[0] = new_nam

        if not a_idx:
            new_nam = 'beta_' + '_'.join(spl)
            string[0] = new_nam

    outStr = ','.join(string)
    out.append(outStr)

with open(str(sys.argv[2]), 'wt') as f:
    f.write("%s\n" % ','.join(fields)+"\n")
    for string in out:
        f.write("%s\n" % string)
