#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
CDR3translator
https://innate2adaptive.github.io/Decombinator/

Take decombined data and translates/extracts the CDR3 sequences.
In order to be classified as (potentially) productive, a rearrangement's CDR3s must be:
    in-frame
    lacking-stop codons
    run from a conserved cysteine to FGXG motif (or appropriate alternatives)

The major change from v3 is that this version exports to the AIRRseq community tsv format, simplifying the process
and crucially giving TCR gene name output in the raw format (in addition to the classic Decombinator fields).

"""

from __future__ import division
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from time import strftime
import argparse
import string
import re
import sys
import collections as coll
import os
import urllib2
import warnings
import gzip

__version__ = '4.1.0'

# Supress Biopython translation warning when translating sequences where length % 3 != 0
warnings.filterwarnings("ignore")

# TODO Potentially add a flag to combine convergent recombinations into a single row?

def args():
    """
    :return: Command line arguments which dictate the script's behaviour
    """

    # Help flag
    parser = argparse.ArgumentParser(
        description='Translate and extract CDR3 sequences from Decombinator classifier files. '
                    'Please see https://innate2adaptive.github.io/Decombinator/ for details.')
    # Add arguments
    parser.add_argument('-in', '--infile', type=str, required=True,
                        help='File containing 5 part Decombinator indexes, (with/without frequencies)')

    parser.add_argument('-c', '--chain', type=str, help='TCR chain (a/b/g/d)', required=False)

    parser.add_argument('-sp', '--species', type=str, required=False, default="human",
                        help='Specify which species TCR repertoire the data consists of (human or mouse). '
                             'Default = human')

    parser.add_argument('-tg', '--tags', type=str, required=False, default="extended",
                        help='Specify which Decombinator tag set to use (extended or original). Default = extended')

    parser.add_argument('-s', '--suppresssummary', action='store_true', required=False,
                        help='Suppress the production of summary data log')

    parser.add_argument('-npf', '--nonproductivefilter', action='store_true', required=False,
                        help='Filter out non-productive reads from the output')

    parser.add_argument('-dz', '--dontgzip', action='store_true', required=False,
                        help='Stop the output FASTQ files automatically being compressed with gzip')

    parser.add_argument('-dc', '--dontcount', action='store_true', required=False,
                        help='Stop printing the running count')

    parser.add_argument('-tfdir', '--tagfastadir', type=str, required=False, default="Decombinator-Tags-FASTAs",
                        help='Path to folder containing TCR FASTA and Decombinator tag files, for offline analysis.'
                             'Default = \"Decombinator-Tags-FASTAs\".')

    return parser.parse_args()


def findfile(filename):
    """
    :param filename: Check whether input file exists or not
    :return: Nothing: script exits if given input file does not exist
    """

    try:
        testopen = open(str(filename), "rU")
        testopen.close()
    except Exception:
        print 'Cannot find the specified input file. Please try again'
        sys.exit()


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
        fl_opener = open

    elif os.path.isfile(expected_dir_name + os.sep + expected_file):
        fl = expected_dir_name + os.sep + expected_file
        fl_opener = open

    else:
        try:
            fl = "https://raw.githubusercontent.com/innate2adaptive/Decombinator-Tags-FASTAs/master/" + expected_file
            urllib2.urlopen(urllib2.Request(fl))  # Request URL, see whether is found
            fl_opener = urllib2.urlopen

        except Exception:
            print "Cannot find following file locally or online:", expected_file
            print "Please either run Decombinator with internet access, or point Decombinator to local copies " \
                  "of the tag and FASTA files with the \'-tf\' flag."
            sys.exit()

    # Return opened file, for either FASTA or tag file parsing
    return fl_opener(fl)


def sort_permissions(fl):
    """
    Need to ensure proper file permissions on output data.
    If users are running pipeline through Docker might otherwise require root access
    :param fl: The file to sort permissions on
    :return: Nothing: script edits permissions where appropriate, if possible
    """

    if oct(os.stat(fl).st_mode)[4:] != '666':
        os.chmod(fl, 0o666)


def import_gene_information(inputargs):
    """
    Obtains gene-specific information for translation
    Runs first: reads in V and J gene sequence and name data (from fasta files)
    and positions of conserved cysteine residues in V genes (from separate files)

    If files cannot be found in local directory, script looks for them online at GitHub

    NB that a number of psuedogenes have no officially designated conserved C (or indeed a 3' C at all)
      Where possible, the nearest suitable C residue is used, where not an arbitrary position of 0 is given
      Moot, as most psuedogenes contain a number of stop codons and thus cannot produce productive rearrangements

    First check that valid tag/species combinations have been used
    :param inputargs: command line (argparse) input arguments dictionary
    :return: Multiple items of TCR data: the V regions (sequence), J regions (sequence), V gene names,
    J gene names, V conserved C translate positions, V conserved position residue identidy, J conserved F
    translate position, J conserved position residue identity, V gene functionality, J gene functionality
    """

    global chainnams, chain
    chain = inputargs['chain']

    if inputargs['tags'] == "extended" and inputargs['species'] == "mouse":
        print "Please note that there is currently no extended tag set for mouse TCR genes.\n" \
              "Decombinator will now switch the tag set in use from \'extended\' to \'original\'.\n" \
              "In future, consider editing the script to change the default, " \
              "or use the appropriate flags (-sp mouse -tg original)."
        inputargs['tags'] = "original"

    if inputargs['tags'] == "extended" and (chain == 'g' or chain == 'd'):
        print "Please note that there is currently no extended tag set for gamma/delta TCR genes.\n" \
              "Decombinator will now switch the tag set in use from \'extended\' to \'original\'.\n" \
              "In future, consider editing the script to change the default, or use the appropriate flags."
        inputargs['tags'] = "original"

    # Check species information
    if inputargs['species'] not in ["human", "mouse"]:
        print "Species not recognised. Please select either \'human\' (default) or \'mouse\'.\n" \
              "If mouse is required by default, consider changing the default value in the script."
        sys.exit()

    # Look for tag and V/J fasta and cysteine position files: if these cannot be found in the working directory,
    # source them from GitHub repositories
    # Note that fasta/tag files fit the pattern "species_tagset_gene.[fasta/tags]"
    # I.e. "[human/mouse]_[extended/original]_TR[A/B/G/D][V/J].[fasta/tags]"

    for gene in ['v', 'j']:
        # Get FASTA data
        fasta_file = read_tcr_file(inputargs['species'], inputargs['tags'], gene, "fasta", inputargs['tagfastadir'])
        globals()[gene + "_genes"] = list(SeqIO.parse(fasta_file, "fasta"))
        fasta_file.close()
        globals()[gene + "_regions"] = [str(string.upper(item.seq)) for item in globals()[gene + "_genes"]]
        globals()[gene + "_names"] = [str(string.upper(item.id).split("|")[1]) for item in globals()[gene + "_genes"]]

        # Get conserved translation residue sites and functionality data
        translation_file = read_tcr_file(inputargs['species'], inputargs['tags'], gene, "translate",
                                         inputargs['tagfastadir'])
        translate_data = [x.rstrip() for x in list(translation_file)]
        translation_file.close()
        globals()[gene + "_translate_position"] = [int(x.split(",")[1]) for x in translate_data]
        globals()[gene + "_translate_residue"] = [x.split(",")[2] for x in translate_data]
        globals()[gene + "_functionality"] = [x.split(",")[3] for x in translate_data]

        if gene == 'v':
            # Get germline CDR data
            cdr_file = read_tcr_file(inputargs['species'], inputargs['tags'], gene, "cdrs", inputargs['tagfastadir'])
            cdr_data = [x.rstrip() for x in list(cdr_file)]
            cdr_file.close()
            v_cdr1 = [x.split(" ")[1] for x in cdr_data]
            v_cdr2 = [x.split(" ")[2] for x in cdr_data]

    return v_regions, j_regions, v_names, j_names, v_translate_position, v_translate_residue, \
           j_translate_position, j_translate_residue, v_functionality, j_functionality, v_cdr1, v_cdr2


def get_cdr3(dcr, headers):
    """
    Checks the productivity of a given DCR-assigned rearrangement.
    Note it requires certain items to be in memory: import_gene_information() must be run first
    :param dcr: the 5 part Decombinator identifier of a given sequence
    :param headers: the headers of the fields that will appear in the final output file (including empty ones)
    :return: a dictionary of the relevant output fields, for downstream transcription into the out file
    """

    # NB: A productively rearranged receptor does not necessarily mean that it is the working receptor used in a cell!
    out_data = coll.defaultdict()
    for field in headers:
        out_data[field] = ''

    out_data['decombinator_id'] = dcr
    out_data['rev_comp'] = 'F'

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

    # TODO remove 'split' if and when the gene names in the tag files get properly adjusted to be consistent
    out_data['v_call'] = v_names[v].split('*')[0]
    out_data['j_call'] = j_names[j].split('*')[0]

    if vdel == 0:
        v_used = v_regions[v]
    else:
        v_used = v_regions[v][:-vdel]

    j_used = j_regions[j][jdel:]

    out_data['sequence'] = ''.join([v_used, ins_nt, j_used])

    # 2. Translate
    out_data['sequence_aa'] = str(Seq(out_data['sequence'], generic_dna).translate())

    # 3. Check whether whole rearrangement is in frame
    if (len(out_data['sequence']) - 1) % 3 == 0:
        out_data['productive'] = 'T'
        out_data['vj_in_frame'] = 'T'
    else:
        out_data['productive'] = 'F'
        out_data['vj_in_frame'] = 'F'

    # 4. Check for stop codons in the in-frame rearrangements
    if '*' in out_data['sequence_aa']:
        out_data['productive'] = 'F'
        out_data['stop_codon'] = 'T'
    else:
        out_data['stop_codon'] = 'F'

    # 5. Check for conserved cysteine in the V gene
    if out_data['sequence_aa'][v_translate_position[v] - 1] == v_translate_residue[v]:
        start_cdr3 = v_translate_position[v] - 1
        out_data['conserved_c'] = 'T'
    else:
        out_data['productive'] = 'F'
        out_data['conserved_c'] = 'F'

    # 5.5 Having found conserved cysteine, only need look downstream to find other end of CDR3
    downstream_c = out_data['sequence_aa'][start_cdr3:]

    # 6. Check for presence of FGXG motif (or equivalent)
    site = downstream_c[j_translate_position[j]:j_translate_position[j] + 4]

    if re.findall(j_translate_residue[j], site):
        end_cdr3 = len(downstream_c) + j_translate_position[j] + start_cdr3 + 1
        out_data['conserved_f'] = 'T'
    else:
        out_data['productive'] = 'F'
        out_data['conserved_f'] = 'F'

    if out_data['productive'] == 'T':
        out_data['junction_aa'] = out_data['sequence_aa'][start_cdr3:end_cdr3]
        out_data['junction'] = out_data['sequence'][start_cdr3 * 3:3 * end_cdr3]
        out_data['cdr1_aa'] = v_cdr1[v]
        out_data['cdr2_aa'] = v_cdr2[v]

    return out_data


out_headers = ['sequence_id', 'v_call', 'd_call', 'j_call', 'junction_aa', 'duplicate_count', 'sequence',
               'junction', 'decombinator_id', 'rev_comp', 'productive', 'sequence_aa', 'cdr1_aa', 'cdr2_aa',
               'vj_in_frame', 'stop_codon', 'conserved_c', 'conserved_f',
               'sequence_alignment', 'germline_alignment', 'v_cigar', 'd_cigar', 'j_cigar']


if __name__ == '__main__':

    # Check input files and parameters
    inputargs = vars(args())
    counts = coll.Counter()

    if inputargs['infile'].endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    # Get chain information
    if not inputargs['chain']:
        # If chain not given, try and infer from input file name
        chaincheck = [x for x in ["alpha", "beta", "gamma", "delta"] if x in inputargs['infile'].lower()]
        if len(chaincheck) == 1:
            chain = chaincheck[0][0]
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

    inputargs['chain'] = chain  # Correct inputarg chain value so that import gene function gets correct input

    suffix = ".tsv"

    filename = inputargs['infile']
    findfile(filename)

    # Extract CDR3s
    v_regions, j_regions, v_names, j_names, v_translate_position, v_translate_residue, j_translate_position, \
    j_translate_residue, v_functionality, j_functionality, v_cdr1, v_cdr2 = import_gene_information(inputargs)

    infile = opener(filename, "rU")

    counts['line_count'] = 0

    # Count non-productive rearrangments
    chainnams = {"a": "alpha", "b": "beta", "g": "gamma", "d": "delta"}

    print "Translating", chainnams[chain], "chain CDR3s from", inputargs['infile']

    outfilename = filename.split(".")[0] + suffix
    with opener(filename, 'rU') as in_file, open(outfilename, 'w') as out_file:
        out_file.write('\t'.join(out_headers) + '\n')

        for line in in_file:

            counts['line_count'] += 1

            comma = [m.start() for m in re.finditer(',', line)]

            if len(comma) == 4:  # pure DCR file, just the five fields (no frequency)
                in_dcr = line.rstrip()
                use_freq = False
                frequency = 1
            elif len(comma) == 5:  # assume that we're working from a .freq file (or equivalent)
                in_dcr = str(line[:comma[4]])
                frequency = int(line[comma[4] + 2:].rstrip())
                use_freq = True
            elif len(comma) > 5:  # assume that it's an n12 file (no frequency)
                print "Incorrect number of comma-delimited fields detected. Please check input and try again."
                sys.exit()
            else:
                print "Based on the number of commas per line, script is unable to determine file type. " \
                      "Please ensure you're inputting a valid file (e.g. .n12 or .freq)."
                sys.exit()

            cdr3_data = get_cdr3(in_dcr, out_headers)
            cdr3_data['sequence_id'] = str(counts['line_count'])

            v = int(line[:comma[0]])
            j = int(line[comma[0] + 2:comma[1]])

            cdr3_data['duplicate_count'] = frequency

            if cdr3_data['productive'] == 'T':
                counts['prod_recomb'] += 1
                productivity = "P"
                out_file.write('\t'.join([str(cdr3_data[x]) for x in out_headers]) + '\n')
            else:
                productivity = "NP"
                counts['NP_count'] += 1
                if not inputargs['nonproductivefilter']:
                    out_file.write('\t'.join([str(cdr3_data[x]) for x in out_headers]) + '\n')

            # Count the number of number of each type of gene functionality (by IMGT definitions, based on prototypic)
            if inputargs['tags'] == 'extended' and inputargs['species'] == 'human':
                counts[productivity + "_" + "V-" + v_functionality[v]] += 1
                counts[productivity + "_" + "J-" + j_functionality[j]] += 1

    print "CDR3 data written to", outfilename

    # Compress output
    if not inputargs['dontgzip']:
        print "Compressing CDR3 output file to", outfilename + ".gz"

        with open(outfilename) as infile, gzip.open(outfilename + '.gz', 'wb') as outfile:
            outfile.writelines(infile)
        os.unlink(outfilename)

        outfilenam = outfilename + ".gz"

    else:
        outfilenam = outfilename

    sort_permissions(outfilenam)

    # Write data to summary file
    if not inputargs['suppresssummary']:

        # Check for directory and make summary file
        if not os.path.exists('Logs'):
            os.makedirs('Logs')
        date = strftime("%Y_%m_%d")

        # Check for existing date-stamped file
        fname = os.path.splitext(os.path.basename(filename))[0]
        summaryname = "Logs/" + date + "_" + fname + "_CDR3_Translation_Summary.csv"
        if not os.path.exists(summaryname):
            summaryfile = open(summaryname, "w")
        else:
            # If one exists, start an incremental day stamp
            for i in range(2, 10000):
                summaryname = "Logs/" + date + "_" + fname + \
                              "_CDR3_Translation_Summary" + str(i) + ".csv"
                if not os.path.exists(summaryname):
                    summaryfile = open(summaryname, "w")
                    break

        # Generate string to write to summary file
        summstr = "Property,Value\nDirectory," + os.getcwd() + "\nInputFile," \
                  + inputargs['infile'] + "\nOutputFile," + outfilenam \
                  + "\nDateFinished," + date + "\nTimeFinished," \
                  + strftime("%H:%M:%S") + "\n\nInputArguments:,\n"
        for s in ['species', 'chain', 'tags', 'dontgzip']:
            summstr = summstr + s + "," + str(inputargs[s]) + "\n"

        summstr = summstr + "\nNumberUniqueDCRsInput," + str(counts['line_count']) \
                  + "\nNumberUniqueDCRsProductive," + str(counts['prod_recomb']) \
                  + "\nNumberUniqueDCRsNonProductive," + str(counts['NP_count'])

        if inputargs['tags'] == 'extended' and inputargs['species'] == 'human':
            summstr = summstr + "\n\nFunctionalityOfGermlineGenesUsed,"
            for p in ['P', 'NP']:
                for g in ['V', 'J']:
                    for f in ['F', 'ORF', 'P']:
                        target = p + '_' + g + '-' + f
                        summstr = summstr + '\n' + target + ',' + str(counts[target])

        print >> summaryfile, summstr
        summaryfile.close()
        sort_permissions(summaryname)

    sys.exit()
