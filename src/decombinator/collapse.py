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

  NOTE - from version 4.1 this optionally looks for barcode 6NI86N at the beginning of the read; instead of M13_6N_I8_6N_I8
  (i.e. only one spacer).
  This makes it compatible with the multiplex protocol in which the barcode is incorproated in the RT step
  In order to work, you must specify an additional command line parameter -ol (see below)
##################
###### INPUT #####
##################

  Required inputs
  -in/--infile : Defines input file.   Takes as input .n12 files produced by Decombinator (v3 or higher),
                    assuming it has been run on suitably barcoded and demultiplexed data.
  -ol/--oligo  : Specifies the spacer (protocol dependent) as M13, I8, I8_single, NEBIO, or TAKARA. The I8 protocol is deprecated.

  Other optional flags:

    -s/--supresssummary: Suppress the production of a summary file containing details of the run into a 'Logs' directory.

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

#######################################################################################################################################################
"""

import ast
import collections as coll
import gzip
import os
import sys
import time
import typing
from importlib import metadata

import networkx as nx
import polyleven
import pyrepseq.nn as prsnn
import regex
import scipy.sparse
from scipy import sparse

########################################################################################################################
# Functions


def num_check(poss_int):
    """Check whether string is feasibly an integer value of zero or greater"""
    try:
        if int(poss_int) >= 0:
            return True
        else:
            return False
    except ValueError:
        return False


def is_dna(poss_dna):
    """Check whether string is feasibly a DNA sequence read"""
    return set(poss_dna.upper()).issubset({"A", "C", "G", "T", "N"})


def check_dcr_file(infile, opener):
    """
    Perform sanity check on input Decombinator file
    Check first few lines to see whether they fit the correct criteria this script relies on
    """
    # Check whether file a) exists and b) is not empty
    if os.path.isfile(infile) == False:
        print("Cannot find file, please double-check path.")
        return False
    print(os.path.getsize(infile))
    if os.path.getsize(infile) == 0:
        raise ValueError(
            "Input file appears to be empty; please double-check path."
        )

    # Check first few lines
    with opener(infile, "rt") as poss_dcr:
        for i in range(5):

            # Check it's a comma-delimited file
            try:
                nextline = next(poss_dcr)
            except StopIteration:
                raise StopIteration(
                    f"Input Decombinator file sanity check warning: {i} line(s) in input file."
                )

            if "," not in nextline:
                print(
                    "Input Decombinator file sanity check fail: seemingly not comma-delimited text file."
                )
                return False
            else:
                fields = nextline.rstrip().split(", ")

            # TODO: reimplement these checks as tests
            # Check lines contain correct number of fields
            if len(fields) != 10:
                print(
                    "Input Decombinator file sanity check fail: file does not contain the correct number of comma-delimited fields (ten)."
                )
                return False
            # Check DCR classifiers are feasible (i.e. 4x integers followed by a DNA string)
            elif not all([num_check(field) for field in fields[0:4]]):
                print(
                    "Input Decombinator file sanity check fail: integer components of Decombinator classifier not feasible (i.e. not integers >= zero)."
                )
                return False
            elif not is_dna(fields[4]):
                print(
                    "Input Decombinator file sanity check fail: Decombinator insert sequence field contains non-DNA sequence."
                )
                return False
            # Check inter-tag and barcode sequences are legitimate DNA strings
            elif is_dna(fields[6]) != True or is_dna(fields[8]) != True:
                print(
                    "Input Decombinator file sanity check fail: inter-tag and/or barcode sequence contains non-DNA sequences."
                )
                return False
            # Check inter-tag and barcode quality are feasible FASTQ quality scores
            elif (
                all(set([num_check(x) for x in get_qual_scores(fields[7])]))
                != True
                or all(set([num_check(x) for x in get_qual_scores(fields[9])]))
                != True
            ):
                print(
                    "Input Decombinator file sanity check fail: inter-tag and/or barcode quality strings do not appear to contain valid scores."
                )
                return False
            # Check that inter-tag and barcode sequence/quality pairs are the correct length
            elif len(fields[6]) != len(fields[7]) or len(fields[8]) != len(
                fields[9]
            ):
                print(
                    "Input Decombinator file sanity check fail: inter-tag and/or barcode sequence and quality string pairs are not of the same length."
                )
                return False

        return True  # If first few lines all pass, assume the file is fine and carry on with the analysis.


def getOligo(oligo_name):
    # New oligos can be added here, specifying their spacers in the given format, and adding them to
    # the returned list.
    oligos = {}
    oligos["m13"] = {"spcr1": "GTCGTGACTGGGAAAACCCTGG", "spcr2": "GTCGTGAT"}
    oligos["i8"] = {"spcr1": "GTCGTGAT", "spcr2": "GTCGTGAT"}
    oligos["i8_single"] = {"spcr1": "ATCACGAC"}
    oligos["nebio"] = {"spcr1": "TACGGG"}
    oligos["takara"] = {"spcr1": "GTACGGG"}

    if oligo_name.lower() not in oligos:
        print(
            "Error: Failed to recognise oligo name. Please choose from "
            + str(list(oligos.keys()))
        )
        sys.exit()

    return oligos[oligo_name.lower()]


def findSubs(subseq, seq):
    # allow up to two substitutions in subseq.
    err_subseqs = regex.findall("(" + subseq + "){1s<=2}", seq)
    return err_subseqs


def findSubsInsOrDels(subseq, seq):
    # allow up to two subsitutions OR (one deletion or one insertion) in subseq.
    err_subseqs = regex.findall("(" + subseq + "){2i+2d+1s<=2}", seq)
    return err_subseqs


def spacerSearch(subseq, seq):
    # first search for subseq within seq. If unsuccessful, allow for subsitutions in subseq.
    # If still unsuccessful, allow for deletions, insertions and subsitutions in subseq.
    foundseq = regex.findall(subseq, seq)
    if not foundseq:
        foundseq = findSubs(subseq, seq)
    if not foundseq:
        foundseq = findSubsInsOrDels(subseq, seq)
    return foundseq


def findFirstSpacer(oligo, seq, oligo_start, oligo_end):
    spacer = []
    spcr1 = oligo["spcr1"]
    spacer += spacerSearch(spcr1, seq[oligo_start:oligo_end])
    return spacer


def findSecondSpacer(oligo, seq):
    spacer = []
    spcr1 = oligo["spcr1"]
    spcr2 = oligo["spcr2"]
    spacer += spacerSearch(spcr2, seq[len(spcr1) :])
    return spacer


def getSpacerPositions(bcseq, spacers):
    # locates the positions of the determined spacers in the sequence
    positions = []
    startpos = 0
    for x in spacers:
        positions.append(bcseq.find(x, startpos))
        startpos += len(x)
    return positions


def filterShortandLongBarcodes(
    b1len: int, b2end: int, bcseq: str, counts: dict
) -> bool:
    if b1len <= 3:
        counts["getbarcode_fail_n1tooshort"] += 1
        return False
    elif b1len >= 9:
        counts["getbarcode_fail_n1toolong"] += 1
        return False
    elif b2end > len(bcseq):
        counts["getbarcode_fail_n2pastend"] += 1
        return False
    else:
        return True


def logExactOrRegexMatch(spacers: list[str], oligo: dict[str], counts: dict):
    if spacers == list(oligo.values()):
        counts["getbarcode_pass_exactmatch"] += 1
    else:
        counts["getbarcode_pass_regexmatch"] += 1


def logFuzzyMatching(
    b1len: int,
    bclength: int,
    spacers: list[str],
    oligo: dict[str],
    counts: dict,
) -> None:
    if b1len == bclength and spacers != list(oligo.values()):
        counts["getbarcode_pass_fuzzymatch_rightlen"] += 1
    elif b1len in [4, 5] and spacers != list(oligo.values()):
        counts["getbarcode_pass_fuzzymatch_short"] += 1
    elif b1len >= 7 and spacers != list(oligo.values()):
        counts["getbarcode_pass_fuzzymatch_long"] += 1
    elif b1len == bclength:
        counts["getbarcode_pass_other"] += 1


def set_barcode(
    fields: list[str], bc_locs: list[int], inputargs: dict
) -> tuple[str, str]:
    # account for N1 barcode being greater or shorter than 6 nt (due to manufacturing errors)
    if str.lower(inputargs["oligo"]) in ["nebio", "takara"]:
        barcode = fields[8][bc_locs[0] : bc_locs[1]]
        barcode_qualstring = fields[9][bc_locs[0] : bc_locs[1]]
    else:
        if (bc_locs[1] - bc_locs[0]) == 6:
            barcode = (
                fields[8][bc_locs[0] : bc_locs[1]]
                + fields[8][bc_locs[2] : bc_locs[3]]
            )
            barcode_qualstring = (
                fields[9][bc_locs[0] : bc_locs[1]]
                + fields[9][bc_locs[2] : bc_locs[3]]
            )

        elif (bc_locs[1] - bc_locs[0]) < 6:
            n1_diff_len = 6 - (bc_locs[1] - bc_locs[0])
            barcode = (
                fields[8][bc_locs[0] : bc_locs[1]]
                + "S" * n1_diff_len
                + fields[8][bc_locs[2] : bc_locs[3]]
            )
            barcode_qualstring = (
                fields[9][bc_locs[0] : bc_locs[1]]
                + "?" * n1_diff_len
                + fields[9][bc_locs[2] : bc_locs[3]]
            )

            counts["readdata_short_barcode"] += 1

        elif (bc_locs[1] - bc_locs[0]) > 6:
            n1_diff_len = 6 - (bc_locs[1] - bc_locs[0])
            barcode = (
                fields[8][bc_locs[0] : bc_locs[0] + 5]
                + "L"
                + fields[8][bc_locs[2] : bc_locs[3]]
            )
            barcode_qualstring = (
                fields[9][bc_locs[0] : bc_locs[0] + 5]
                + "?" * n1_diff_len
                + fields[9][bc_locs[2] : bc_locs[3]]
            )

            counts["readdata_long_barcode"] += 1

    # L and S characters get quality scores of "?", representative of Q30 scores
    return barcode, barcode_qualstring


def get_qual_scores(qualstring):
    """Convert FASTQ quality scores to their integer Q score equivalent"""
    quality_list = [ord(x) - 33 for x in qualstring]
    return quality_list


def get_err_prob(Q):
    """Returns the probability of a given base call being incorrect, based on quality score"""
    return 10 ** (-Q / 10)


def check_umi_quality(qualstring: tuple[str], parameters: list[int]) -> bool:
    """
    Input: string representing barcode quality and quality check parameters
    Output: False if barcode fails quality check
            True if barcode passes quality check
    """
    quality_list = get_qual_scores(qualstring)
    number_below_min = sum([x < parameters[0] for x in quality_list])
    average_quality = sum(quality_list) / len(quality_list)
    return number_below_min > parameters[1] or average_quality < parameters[2]


def are_seqs_equivalent(seq1, seq2, lev_percent_threshold):
    # Returns True if seqs can be considered the same, False otherwise
    # Definition of equivalent:
    #   levenshtein distance as a percentage of the shorter of the two seqs is <= threshold
    threshold = len(min(seq1, seq2, key=len)) * lev_percent_threshold
    return polyleven.levenshtein(seq1, seq2) <= threshold


def are_barcodes_equivalent(bc1, bc2, threshold):
    return polyleven.levenshtein(bc1, bc2) <= threshold


def get_barcode_positions(
    bcseq: str, inputargs: dict, counts: dict
) -> list[int]:
    """
    Given a barcode-region sequence, outputs the sequence of the do-docamer barcode.
    For m13 and i8 oligos, this barcode (theoretically) consists of the concatentation of the two random hexamer sequences contained in the ligation oligo.
    However errors in sequences and ligation oligo production can mean that the random nucleotides are not always at the expected position.
    This function uses the known sequence of the spacers (which bound each of the two N6s to their 5') to deduce the random sequences.
    Returns a list of four numbers, giving the start and stop positions of N1 and N2 respectively.
    """
    if str.lower(inputargs["oligo"]) not in [
        "i8",
        "i8_single",
        "m13",
        "nebio",
        "takara",
    ]:
        raise ValueError(
            "The flag for the -ol input must be one of M13, I8, I8_single, NEBIO, or TAKARA."
        )

    if (
        "N" in bcseq and inputargs["allowNs"] == False
    ):  # ambiguous base-call check
        counts["getbarcode_fail_N"] += 1
        return

    # gets spacer sequences of specified oligo
    oligo = getOligo(inputargs["oligo"])

    # sets first spacer based on specified oligo
    if str.lower(inputargs["oligo"]) == "nebio":
        oligo_start = 18
        oligo_end = oligo_start + 10
    elif str.lower(inputargs["oligo"]) == "takara":
        oligo_start = 0
        oligo_end = 19
    else:
        oligo_start = 0
        allowance = 10
        oligo_end = allowance + len(oligo["spcr1"])
    spacers = findFirstSpacer(oligo, bcseq, oligo_start, oligo_end)

    # sequences with no first spacer are removed from analysis
    if not len(spacers) == 1:
        counts["getbarcode_fail_nospacerfound"] += 1
        return None

    # sets second spacer based on specified oligo (unless single oligo)

    if str.lower(inputargs["oligo"]) not in ["i8_single", "nebio", "takara"]:
        spacers += findSecondSpacer(oligo, bcseq)
        # sequences which do not have two spacers are logged then removed from analysis
        if not len(spacers) == 2:
            counts["getbarcode_fail_not2spacersfound"] += 1
            return None
    spacer_positions = getSpacerPositions(bcseq, spacers)

    if str.lower(inputargs["oligo"]) in ["nebio", "takara"]:
        # set expected barcode length
        if str.lower(inputargs["oligo"]) == "nebio":
            bclength = 17
        else:
            bclength = 12
        # start and end of barcode positions are set
        b1start = 0
        b1end = bclength

        b1len = b1end - b1start

        # filtering and logging
        logExactOrRegexMatch(spacers, oligo, counts)
        logFuzzyMatching(b1len, bclength, spacers, oligo, counts)

        return [b1start, b1end]

    else:
        if str.lower(inputargs["oligo"]) == "i8_single":
            # set expected barcode length
            bclength = 6
            # start and end of barcode positions are set
            b1start = 0
            b1end = spacer_positions[0]
            b2start = spacer_positions[0] + len(spacers[0])

            b2end = b2start + bclength
            b1len = b1end - b1start

            # filtering and logging
            if not filterShortandLongBarcodes(b1len, b2end, bcseq, counts):
                return None
            logExactOrRegexMatch(spacers, oligo, counts)
            logFuzzyMatching(b1len, bclength, spacers, oligo, counts)

            return [b1start, b1end, b2start, b2end]

        else:
            # set expected barcode length
            bclength = 6
            # start and end of barcode positions are set
            b1start = spacer_positions[0] + len(spacers[0])
            b1end = spacer_positions[1]
            b2start = spacer_positions[1] + len(spacers[1])
            b2end = b2start + bclength
            b1len = b1end - b1start

            # filtering and logging
            if not filterShortandLongBarcodes(b1len, b2end, bcseq, counts):
                return None
            logExactOrRegexMatch(spacers, oligo, counts)
            logFuzzyMatching(b1len, bclength, spacers, oligo, counts)

            return [b1start, b1end, b2start, b2end]


def read_in_data(
    data,
    inputargs,
    barcode_quality_parameters,
    lev_threshold,
    dont_count,
    opener,
):
    ###########################################
    ############# READING DATA IN #############
    ###########################################

    # Check whether file appears to contain suitable verbose Decombinator output for collapsing
    # TODO: reimplement this section as tests
    if inputargs["command"] == "collapse":
        if not inputargs["dontcheckinput"]:
            if not check_dcr_file(data, opener):
                print(
                    "Please check that file contains suitable Decombinator output for collapsing."
                )
                print(
                    "Alternatively, disable the input file sanity check by changing the 'dontcheckinput' flag, i.e. '-di True'"
                )
                sys.exit()
        data = opener(data, "rt")

    if not data:
        raise ValueError(
            "No reads found in input file. Check .n12 and log files for errors."
        )

    print("Reading data in...")
    t0 = time.time()
    barcode_dcretc = coll.defaultdict(list)
    barcode_lookup = coll.defaultdict(list)

    input_dcr_counts = coll.Counter()
    ratio = 1
    l = 0

    for lcount, line in enumerate(data):
        if inputargs["command"] == "collapse":
            line = line.rstrip("\n").split(", ")
        if ratio < 0.01 and (time.time() - t0) > 3600:
            break
        if lcount % 50000 == 0 and lcount != 0 and not dont_count:
            print(
                "   Read in",
                lcount,
                "lines... ",
                round(time.time() - t0, 2),
                "seconds",
            )
            ratio = (len(barcode_lookup) - l) / len(barcode_lookup)
            l = len(barcode_lookup)
            print(round(ratio, 2))

        counts["readdata_input_dcrs"] += 1

        bc_locs = get_barcode_positions(line[8], inputargs, counts)
        if not bc_locs:
            counts["readdata_fail_no_bclocs"] += 1
            continue

        barcode, barcode_qualstring = set_barcode(line, bc_locs, inputargs)
        # L and S characters get quality scores of "?", representative of Q30 scores

        if check_umi_quality(barcode_qualstring, barcode_quality_parameters):
            # barcode is not sufficient quality, skip to next line of file
            counts["readdata_fail_low_barcode_quality"] += 1
            continue

        dcr = line[:5]
        input_dcr_counts[str(dcr)] += 1

        seq = line[6]

        if len(seq) > inputargs["lenthreshold"]:
            # end V tag to start J tag too long to be sane
            counts["readdata_fail_overlong_intertag_seq"] += 1
            continue

        counts["readdata_success"] += 1
        seq_qualstring = line[7]
        seq_id = line[5]
        dcretc = "|".join([str(dcr), seq, seq_qualstring, seq_id])

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
                    barcode_dcretc[
                        barcode + "|" + str(index[0]) + "|" + index[1]
                    ].append(dcretc)
                    protodcretc_list = barcode_dcretc[
                        barcode + "|" + str(index[0]) + "|" + index[1]
                    ]
                    seq_counter = coll.Counter(
                        map(lambda x: x.split("|")[1], protodcretc_list)
                    )
                    protoseq = seq_counter.most_common(1)[0][
                        0
                    ]  # find most common sequence in group

                    if not index[1] == protoseq:
                        # if there is a new protoseq, replace record with old protoseq
                        # with identical record with updated  protoseq
                        barcode_dcretc[
                            barcode + "|" + str(index[0]) + "|" + protoseq
                        ] = barcode_dcretc[
                            barcode + "|" + str(index[0]) + "|" + index[1]
                        ]
                        del barcode_dcretc[
                            barcode + "|" + str(index[0]) + "|" + index[1]
                        ]

                        barcode_lookup[barcode][index[0]] = [
                            index[0],
                            protoseq,
                        ]

                group_assigned = True
                # if assigned to a group, stop and move onto next read
                break

            if not group_assigned:
                # if no appropriate group found, create new group with correctly incremented index
                barcode_lookup[barcode].append([index[0] + 1, seq])
                barcode_dcretc[
                    "|".join([barcode, str(index[0] + 1), seq])
                ].append(dcretc)
                group_assigned = True

        else:
            # if no identical barcode found, create new barcode group with index zero
            barcode_lookup[barcode].append([0, seq])
            barcode_dcretc["|".join([barcode, "0", seq])].append(dcretc)
            group_assigned = True

    counts["readdata_barcode_dcretc_keys"] = len(barcode_dcretc.keys())
    counts["number_input_unique_dcrs"] = len(input_dcr_counts.keys())
    counts["number_input_total_dcrs"] = sum(input_dcr_counts.values())

    t1 = time.time()
    print("   Read in total of", lcount + 1, "lines")
    print(
        "  ",
        counts["readdata_success"],
        "reads sorted into",
        len(barcode_dcretc),
        "initial groups",
    )
    print("  ", round(t1 - t0, 2), "seconds")
    counts["time_readdata_s"] = t1

    return barcode_dcretc


def create_clustering_objs(
    barcode_dcretc: dict[str, list[str]],
) -> tuple[int, list[tuple[str, str]], list[tuple[str, str]]]:

    # get number of initial groups
    num_initial_groups = len(barcode_dcretc)

    # convert barcode_dcretc collection to list format
    barcode_dcretc_list = []
    for _, (j, k) in enumerate(barcode_dcretc.items()):
        barcode_dcretc_list.append((j, k))

    umi_protoseq_tuple = [
        (x[0].split("|")[0], x[0].split("|")[2]) for x in barcode_dcretc_list
    ]

    return num_initial_groups, barcode_dcretc_list, umi_protoseq_tuple


def make_merge_groups(
    umi_protoseq_tuple: list[tuple[str, str]],
    barcode_threshold: int,
    dont_count: bool,
) -> sparse.coo_matrix:
    # cluster similar UMIs
    umi_list = [x[0] for x in umi_protoseq_tuple]
    if len(umi_list) == 0:
        raise ValueError("No UMIs to cluster, check .n12 file for errors")

    print("Clustering UMIs...")
    print("  ", len(umi_list), "unique UMIs")
    matches = prsnn.symdel(
        umi_list,
        max_edits=barcode_threshold,
        progress=not dont_count,
        output_type="coo_matrix",
    )
    matches = sparse.triu(matches)  # Remove duplicates
    matches.sum_duplicates()  # Efficient method to sort made safe by triu

    print(
        "  ",
        matches.getnnz(),
        "UMIs within edit distance of",
        barcode_threshold,
    )

    return matches


def make_clusters(
    merge_groups: sparse.coo_matrix,
    barcode_dcretc: list[tuple[str, list[str]]],
    seq_threshold: int,
) -> coll.defaultdict[str, list[str]]:
    # Considers clusters as an undirected graph composed of disconnected subgraphs.
    # The nodes of the graph are the initial groups of barcode/protosequences. Edges between nodes
    # describe which inital groups form clusters and should be merged.
    # We form a graph of only those initial groups that feature in merge_groups (i.e. that should be
    # merged with one or more other initial groups). The remaining groups that do not need merging are
    # added to clusters after graph analysis at the end of this function.

    # initialise empty collection
    clusters = coll.defaultdict(list)
    n_merged_UMIs = 0

    # initialise empty graph
    G = nx.Graph()

    percent_seq_threshold = seq_threshold / 100.0

    for i, j in zip(merge_groups.row, merge_groups.col):
        protoseqs = [
            barcode_dcretc[i][0].split("|")[2],
            barcode_dcretc[j][0].split("|")[2],
        ]
        if are_seqs_equivalent(
            protoseqs[0], protoseqs[1], percent_seq_threshold
        ):
            G.add_edge(i, j)
            n_merged_UMIs += 1

    print("    ", n_merged_UMIs, "merged UMIs")

    # extracts subgraphs (clusters) from the full graph
    con_comp = nx.connected_components(G)

    for subgraph in con_comp:
        # get full barcode barcode information of the first node in the subgraph from barcode_dcretc
        # this will be serve as the dictionary key for the cluster
        base_node_barcode = barcode_dcretc[list(subgraph)[0]][0]

        # get the full sequence information of each node in the subgraph from barcode_dcretc and
        # add them to cluster collection with a cluster representative barcode (base_node_barcode)
        for k in list(subgraph):
            clusters[base_node_barcode] += barcode_dcretc[k][1]

    # add remaining barcode/protoseqs that do not need merging to the clusters
    for i, bdcretc in enumerate(barcode_dcretc):
        # if already accounted for in the merged_groups then skip over
        if i in G.nodes:
            continue
        else:
            base_node_barcode = bdcretc[0]
            clusters[base_node_barcode] = bdcretc[1]

    return clusters


def write_clusters(clusters):
    # create directory to store cluster data without overwriting exiting directories
    dirname = "clusters"
    count = 1
    while os.path.isdir(dirname):
        dirname = "clusters" + str(count)
        count += 1
    os.mkdir(dirname)

    print(
        "   Writing clusters to directory: ", os.path.abspath(dirname), "..."
    )
    # write data of each cluster to a separate file and store in clusters directory
    for k in clusters:
        with open(
            dirname + os.sep + "|".join(k.split("|")[:2]) + ".txt", "w"
        ) as ofile:
            for j in clusters[k]:
                print(j, file=ofile)
    return 1


def cluster_UMIs(
    barcode_dcretc: coll.defaultdict[str, list[str]],
    inputargs: dict[str, typing.Union[str, bool, int]],
    barcode_threshold: int,
    seq_threshold: int,
    dont_count: bool,
) -> coll.defaultdict[str, list[str]]:
    # input data of form: {'barcode1|index|protoseq': [dcretc1, dcretc2,...], 'barcode2|index|protoseq|: [dcretc1, dcretc2,...], ...}
    # (see read_in_data function for details)
    # This function merges groups that have both equivalent barcodes and equivalent protoseqs
    # output data of form: {'barcode1|index|protoseq': [dcretc1, dcretc2,...], 'barcode2|index|protoseq|: [dcretc1, dcretc2,...], ...}
    # Output format is same as input format, but after clustering. The protoseq is recalculated when new
    # dcretcs are added to a cluster

    print("Clustering barcode groups...")
    t0 = time.time()

    num_initial_groups, barcode_dcretc_list, umi_protoseq_tuple = (
        create_clustering_objs(barcode_dcretc)
    )

    t0 = time.time()
    matches = make_merge_groups(
        umi_protoseq_tuple, barcode_threshold, dont_count
    )

    print("  ", "comparing TCR sequences of similar UMIs...")
    clusters = make_clusters(matches, barcode_dcretc_list, seq_threshold)

    print(
        "  ",
        num_initial_groups,
        "groups merged into",
        len(clusters),
        "clusters",
    )
    t1 = time.time()
    print("  ", round(t1 - t0, 10), "seconds")

    # dump clusters to separate files if desired
    if inputargs["writeclusters"]:
        write_clusters(clusters)

    return clusters


def collapsinate(
    data,
    inputargs,
    barcode_quality_parameters,
    lev_threshold,
    barcode_distance_threshold,
    outpath,
    file_id,
    dont_count,
    opener=None,
):

    # read in, structure, and quality check input data
    barcode_dcretc = read_in_data(
        data,
        inputargs,
        barcode_quality_parameters,
        lev_threshold,
        dont_count,
        opener,
    )

    # cluster similar UMIs
    clusters = cluster_UMIs(
        barcode_dcretc,
        inputargs,
        barcode_distance_threshold,
        lev_threshold,
        dont_count,
    )

    # collapse (count) UMIs in each cluster and print to output file
    print("Collapsing clusters...")
    t0 = time.time()

    collapsed = coll.Counter()
    cluster_sizes = coll.defaultdict(list)

    for c in clusters:
        protodcr = coll.Counter(
            map(lambda x: x.split("|")[0], clusters[c])
        ).most_common(1)[0][
            0
        ]  # find most common dcr in each cluster
        collapsed[protodcr] += 1
        cluster_sizes[protodcr].append(len(clusters[c]))

    counts["number_output_unique_dcrs"] = len(collapsed)
    counts["number_output_total_dcrs"] = sum(collapsed.values())

    t1 = time.time()
    print("  ", round(t1 - t0, 2), "seconds")

    print("Writing to variable...")
    out_data = []

    average_cluster_size_counter = coll.Counter()
    for dcr, dcr_count in collapsed.items():
        av_clus_size = round(sum(cluster_sizes[dcr]) / dcr_count)
        average_cluster_size_counter[av_clus_size] += 1
        list_dcr = ast.literal_eval(
            dcr
        )  # TODO: keep data in object form throughout
        list_dcr.extend([dcr_count, av_clus_size])
        out_data.append(list_dcr)

    # only need to run this bit if interested in the number of times each barcode is repeated in the data
    if inputargs["barcodeduplication"] == True:
        outfile = outpath + file_id + "_barcode_duplication.txt"
        outhandle = open(outfile, "w")
        for bc, copies in clusters.items():
            barcode_index = "|".join(bc.split("|")[:2])
            print(",".join([barcode_index, str(len(copies))]), file=outhandle)
        outhandle.close()
        print("barcode duplication data saved to", outfile)

    counts["outfilenam"] = "Saved to variable"

    return out_data, collapsed, average_cluster_size_counter


def collapsinator(inputargs: dict, data: list = None) -> list:
    """Function wrapper for Collapsinator"""

    print("Running Collapsinator version", metadata.version("decombinator"))
    if inputargs["extension"] == "n12":
        inputargs["extension"] = "freq"
    suffix = "." + inputargs["extension"]

    if inputargs["infile"].endswith(".gz"):
        opener = gzip.open
    else:
        opener = open

    global counts
    counts = coll.Counter()
    counts["start_time"] = time.time()

    ## this is [min_barcode_nt_quality, max_bc_nts_with_min_quality, min_avg_bc_quality]
    barcode_quality_parameters = [
        inputargs["minbcQ"],
        inputargs["bcQbelowmin"],
        inputargs["avgQthreshold"],
    ]

    ## this is the percentage lev distance that is allowed to determine whether two sequences are equivalent
    lev_threshold = inputargs["percentlevdist"]

    ## this is the number of barcode edits that are allowed to call two barcodes equivalent
    barcode_distance_threshold = inputargs["bcthreshold"]

    if inputargs["command"] == "collapse":
        data = inputargs["infile"]
    outpath = ""

    file_id = inputargs["infile"].split("/")[-1].split(".")[0]

    ## this is a boolean for printing progress of the run to the terminal (False for printing, True for not printing; default = False)
    dont_count = inputargs["dontcount"]

    ########################################

    out_data, collapsed, average_cluster_size_counter = collapsinate(
        data,
        inputargs,
        barcode_quality_parameters,
        lev_threshold,
        barcode_distance_threshold,
        outpath,
        file_id,
        dont_count,
        opener,
    )

    counts["end_time"] = time.time()
    counts["time_taken_total_s"] = counts["end_time"] - counts["start_time"]

    #######################################

    # Write data to summary file
    chain = inputargs["chain"]
    chainnams = {"a": "alpha", "b": "beta", "g": "gamma", "d": "delta"}
    if inputargs["suppresssummary"] == False:

        logpath = inputargs["outpath"] + f"Logs{os.sep}"
        file_id_split = file_id.split(os.sep)
        sample_name = file_id_split[-1]

        # Check for directory and make summary file
        if not os.path.exists(logpath):
            os.makedirs(logpath)
        date = time.strftime("%Y_%m_%d")

        # Check for existing date-stamped file
        summaryname = (
            logpath
            + date
            + "_"
            + "dcr_"
            + sample_name
            + f"_{chainnams[chain.lower()]}"
            + "_Collapsing_Summary.csv"
        )
        if not os.path.exists(summaryname):
            summaryfile = open(summaryname, "w")
        else:
            # If one exists, start an incremental day stamp
            for i in range(2, 10000):
                summaryname = (
                    logpath
                    + date
                    + "_"
                    + "dcr_"
                    + sample_name
                    + f"_{chainnams[chain.lower()]}"
                    + "_Collapsing_Summary"
                    + str(i)
                    + ".csv"
                )
                if not os.path.exists(summaryname):
                    summaryfile = open(summaryname, "w")
                    break

        inout_name = (
            "_".join(f"{file_id}".split("_")[:-1])
            + f"_{chainnams[chain.lower()]}"
        )

        # Generate string to write to summary file
        summstr = (
            "Property,Value\nVersion,"
            + str(metadata.version("decombinator"))
            + "\nDirectory,"
            + os.getcwd()
            + "\nInputFile,"
            + inout_name
            + "\nOutputFile,"
            + inout_name
            + "\nDateFinished,"
            + date
            + "\nTimeFinished,"
            + time.strftime("%H:%M:%S")
            + "\nTimeTaken(Seconds),"
            + str(round(counts["time_taken_total_s"], 2))
            + "\n\n"
        )

        for s in [
            "extension",
            "dontgzip",
            "allowNs",
            "dontcheckinput",
            "barcodeduplication",
            "minbcQ",
            "bcQbelowmin",
            "bcthreshold",
            "lenthreshold",
            "percentlevdist",
            "avgQthreshold",
            "positionalbarcodes",
            "oligo",
        ]:
            summstr = summstr + s + "," + str(inputargs[s]) + "\n"

        counts["pc_input_dcrs"] = (
            counts["number_input_total_dcrs"] / counts["readdata_input_dcrs"]
        )
        counts["pc_uniq_dcr_kept"] = (
            counts["number_output_unique_dcrs"]
            / counts["number_input_unique_dcrs"]
        )
        counts["pc_total_dcr_kept"] = (
            counts["number_output_total_dcrs"]
            / counts["number_input_total_dcrs"]
        )

        counts["avg_input_tcr_size"] = (
            counts["number_input_total_dcrs"]
            / counts["number_input_unique_dcrs"]
        )
        counts["avg_output_tcr_size"] = (
            counts["number_output_total_dcrs"]
            / counts["number_output_unique_dcrs"]
        )
        counts["avg_RNA_duplication"] = 1 / counts["pc_total_dcr_kept"]

        # success properties
        summstr = (
            summstr
            + "\nInputUncollapsedDCRLines,"
            + str(counts["readdata_input_dcrs"])
            + "\nUniqueDCRsPassingFilters,"
            + str(counts["number_input_unique_dcrs"])
            + "\nTotalDCRsPassingFilters,"
            + str(counts["number_input_total_dcrs"])
            + "\nPercentDCRPassingFilters(withbarcode),"
            + str(round(counts["pc_input_dcrs"], 3))
            + "\nUniqueDCRsPostCollapsing,"
            + str(counts["number_output_unique_dcrs"])
            + "\nTotalDCRsPostCollapsing,"
            + str(counts["number_output_total_dcrs"])
            + "\nPercentUniqueDCRsKept,"
            + str(round(counts["pc_uniq_dcr_kept"], 3))
            + "\nPercentTotalDCRsKept,"
            + str(round(counts["pc_total_dcr_kept"], 3))
            + "\nAverageInputTCRAbundance,"
            + str(round(counts["avg_input_tcr_size"], 3))
            + "\nAverageOutputTCRAbundance,"
            + str(round(counts["avg_output_tcr_size"], 3))
            + "\nAverageRNAduplication,"
            + str(round(counts["avg_RNA_duplication"], 3))
            + "\n\nBarcodeFail_ContainedNs,"
            + str(counts["getbarcode_fail_N"])
            + "\nBarcodeFail_SpacersNotFound,"
            + str(counts["readdata_fail_no_bclocs"])
            + "\nBarcodeFail_LowQuality,"
            + str(counts["readdata_fail_low_barcode_quality"])
        )

        print(summstr, file=summaryfile)
        summaryfile.close()

        # create output data for UMI histogram and save to file (optional)
        if inputargs["UMIhistogram"]:

            hfileprefix = "_".join(
                summaryname.split("_")[:-2] + ["UMIhistogram"]
            )
            # iterate to create unique file name
            if os.path.exists(hfileprefix + ".csv"):
                i = 1
                while os.path.exists(hfileprefix + str(i) + ".csv"):
                    i += 1
                hfileprefix += str(i)

            hfilename = hfileprefix + ".csv"

            with open(hfilename, "w") as hfile:

                for av, count in sorted(average_cluster_size_counter.items()):
                    print(str(av) + "," + str(count), file=hfile)

            # print instructions for creating histogram plot using script in Supplementary-Scripts
            print(
                "\nAverage UMI cluster size histogram data saved to", hfilename
            )
            print(
                "To plot histogram, please use UMIhistogram.py script located in the Decombinator-Tools repository."
            )
            print(
                "Decombinator-Tools can be found at https://github.com/innate2adaptive/Decombinator-Tools."
            )
            print("With the Decombinator-Tools repository downloaded, run:")
            codestr = (
                "python path/to/Decombinator-Tools/UMIHistogram.py -in "
                + hfilename
            )
            print("#" * (len(codestr) + 4))
            print(" ", codestr, " ")
            print("#" * (len(codestr) + 4))

    del counts
    return out_data


if __name__ == "__main__":
    print(
        "Calling Collapsinator from the shell has been depreciated as of Decombinator V4.3. \
        Please check the README on how to update your script, or alternativley change branch to decombinator_v4.2 \
        which retains this functionality."
    )
