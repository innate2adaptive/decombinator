# Decombinator 4.01
# James M. Heather, August 2016, UCL
# https://innate2adaptive.github.io/Decombinator/

##################
### BACKGROUND ###
##################

# Searches FASTQ reads (produced through Demultiplexor.py) for rearranged TCR chains
# Can currently analyse human and mouse TCRs, both alpha/beta and gamma/delta chains
# NB: Human alpha/beta TCRs are the most thoroughly tested, due to the nature of the data we generated. YMMV.

# The current version (v4) is optimised for interpretation of data generated using our wet lab protocol, but could be modified to work on any data.

# This version incorporates minimal changes and fixes from Decombinator v3.1, which in turn represented improvements upon a combination of the two previously in use Decombinator versions
# i.e. Decombinator V2.2 (written primarily by Nic Thomas, see Thomas et al, Bioinformatics 2013, DOI: 10.1093/bioinformatics/btt004)
# and vDCR (which was v1.4 modified by James Heather, see Heather et al, Frontiers in Immunology 2016, DOI: 10.3389/fimmu.2015.00644)
# Now faster, more accurate and easier to use than either of the previous versions.

# Version 4.2 introduces some changes to work with Demultiplexor 4.01 and to work for the new multiplex protocol. It can look for a
# barcode in either read 1 (multiplex protocol) or read 2 (ligation protocol). This is controled by a new required flag -bc_read
# which must be R1 or R2. The bc_length can also be set - the default is 42

##################
###### INPUT #####
##################

# As with entire pipeline, Decombintator is run using command line arguments to provide user parameters
# All arguments can be read by viewing the help data, by running python Decombintator.py -h

# The two required parameters are
# 1. -fq/--fastq which identify FASTQ reads produced by Demultiplexor.py (unzipped or gzipped).
# 2 -br/--bc_read which determines whether the barcode is obtained from the beginning of R2 as in the standard
# 5'RACE ligation protocol, using the M13-I8-6N-I8-6N or the older SP2-I8-6N-I8-6N oligonucleotide; or from the beginning of R1
# as in the new multiplex protocol. The lenght of sequence conttaining the barcode can also be determeind using -bl/--bc_length.
# The default is 42; but for the multiplex, 22 is enough.

# The TCR chain locus to look for can be explicitly specified using the -c flag
# Users can use their choice of chain identifiers from this list (case insensitive): a/b/g/d/alpha/beta/gamma/delta/TRA/TRB/TRG/TRD/TCRA/TCRB/TCRG/TCRD
# If no chain is provided (or if users which to minimise input arguments), script can infer chain from the FASTQ filename
# I.e. "alpha_sample.fq" would be searched for alpha chain recombinations
# NB: This autodetection only works if there is only ONE TCR locus present in name (which must be spelt out in full)

# Other optional flags:

# -s/--supresssummary: Suppress the production of a summary file containing details of the run into a 'Logs' directory.

# -dz/--dontgzip: Suppress the automatic compression of output demultiplexed FASTQ files with gzip.
# Using this flag makes the script execute faster, but data will require more storage space.

# -dc/--dontcount: Suppress whether or not to show the running line count, every 100,000 reads.
# Helps in monitoring progress of large batches.

# -dk/--dontcheck: Suppress the FASTQ sanity check.
# Strongly recommended to leave alone: sanity check inspects first FASTQ read for basic FASTQ parameters.

# -pf/--prefix: Allows users to specify the prefix of the Decombinator TCR index files produced. Default = 'dcr_'

# -ex/--extension: Allows users to specify the file extension of the Decombinator TCR index files produced. Default = '.n12'

# -or/--orientation: Allows users to specify which DNA orientations to check for TCR reads. Default = reverse only, as that's what the protocol produces.
# This will likely need to be changed for analysing data produced by protocols other than our own.

# -tg/--tags: Allows users to specify which tag set they wish to use. For human alpha/beta TCRs, a new 'extended' tag set is recommended, as it covers more genes.
# Unfortunately an extended tag set is only currently available for human a/b genes.

# -sp/--species: Current options are only human or mouse. Help could potentially be provided for generation of tags for different species upon request.

# -N/--allowNs: Provides users the option to allow 'N's (ambiguous base calls), overriding the filter that typically removes rearrangements that contain them.
# Users are recommended to not allow Ns, as such bases are both themselves low quality data and predict reads that are generally less trustworthy.

# -ln/--lenthreshold: The length threshold which (the inter-tag region of) successful rearrangements must be under to be accepted. Default = 130.

# -tfdir/--tagfastadir: The path to a local copy of a folder containing the FASTA and Decombinator tag files required for offline analysis.
# Ordinarily such files can be downloaded on the fly, reducing local clutter.
# By default the script looks for the required files in the present working directory, then in a subdirectory called "Decombinator-Tags-FASTAs", then online.
# Files are hosted on GitHub, here: https://github.com/innate2adaptive/Decombinator-Tags-FASTAs

# -nbc/--nobarcoding: Run Decombinator without any barcoding, i.e. use the whole read.
# Recommended when running on data not produced using the Innate2Adaptive lab's ligation-mediated amplification protocol

# -bl/--bclength : sets the length of sequence to be stored by Decombinator from R1 or R2 (as set by -bc_read) for further use by Collapsinator.


##################
##### OUTPUT #####
##################

# Produces a '.n12' file by default, which is a standard comma-delimited Decombinator output file with several additional fields:
# V index, J index, # V deletions, # J deletions, insert, ID, TCR sequence, TCR quality, barcode sequence, barcode quality
# NB The TCR sequence given here is the 'inter-tag' region, i.e. the sequence between the start of the found V tag the end of the found J tag

##################
#### PACKAGES ####
##################
from __future__ import division

import collections as coll
import gzip
import itertools
import os
import string
import sys
import urllib
from importlib import metadata
from time import strftime, time

import Levenshtein as lev
from acora import AcoraBuilder
from Bio import SeqIO
from Bio.Seq import Seq

##########################################################
############# FASTQ SANITY CHECK AND PARSING #############
##########################################################


def opener_check(inputargs):
    # Determine compression status (and thus opener required)
    if inputargs["infile"].endswith(".gz"):
        return gzip.open
    else:
        return open


def fastq_check(inputargs, opener, samplenam, summaryname, logpath) -> None:
    """fastq_check(file): Performs a rudimentary sanity check to see whether a file is indeed a FASTQ file"""

    with opener(inputargs["infile"], "rt") as possfq:
        # islice used to very quickly check files
        if sum(1 for _ in itertools.islice(possfq, 4)) < 4:
            # Log of empty file required for pipeline
            if inputargs["suppresssummary"] == False:
                inout_name = (
                    "_".join(f"{samplenam}".split("_")[:-1])
                    + f"_{chainnams[chain]}"
                )
                summstr = (
                    "OutputFile," + inout_name + "\nNumberReadsInput," + "0"
                )
            if not os.path.exists(summaryname):
                summaryfile = open(summaryname, "wt")
            else:
                # If one exists, start an incremental day stamp
                date = strftime("%Y_%m_%d")
                for i in range(2, 10000):
                    summaryname = logpath + date + "_"
                    if inputargs["chain"]:
                        summaryname += chainnams[chain] + "_"
                    summaryname += (
                        samplenam + "_Decombinator_Summary" + str(i) + ".csv"
                    )
                    if not os.path.exists(summaryname):
                        summaryfile = open(summaryname, "wt")
                        break
            print(summstr, file=summaryfile)
            summaryfile.close()
            sort_permissions(summaryname)
            raise ValueError(
                "There are fewer than four lines in this file, and thus it is not a valid FASTQ file. Please check input and try again."
            )
        else:
            read = [i for i in itertools.islice(possfq, 0, 4)]

    # @ check
    if read[0][0] != "@":
        raise ValueError(
            f"Expected @ symbol at beginning of file for valid FASTQ. Found {read[0][0]}."
        )
    # Descriptor check
    if read[2][0] != "+":
        raise ValueError(
            f"Expected + symbol at beginning of third line for valid FASTQ. Found {read[2][0]}."
        )
    # Read/quality match check
    if len(read[1]) != len(read[3]):
        raise ValueError(
            f"Length of read to match length of read quality. Found read length = {len(read[1])} and read quality length = {len(read[3])}"
        )


def revcomp(read):
    """rc(read): Wrapper for SeqIO reverse complement function"""
    return str(Seq(read).reverse_complement())


def read_tcr_file(species, tagset, gene, filetype, expected_dir_name):
    """Reads in the FASTA and tag data for the appropriate TCR locus"""

    # Define expected file name
    expected_file = (
        species
        + "_"
        + tagset
        + "_"
        + "TR"
        + chain.upper()
        + gene.upper()
        + "."
        + filetype
    )

    # First check whether the files are available locally (in pwd or in bundled directory)
    if os.path.isfile(expected_file):
        fl = expected_file
    elif os.path.isfile(expected_dir_name + os.sep + expected_file):
        fl = expected_dir_name + os.sep + expected_file
    else:
        try:
            fl = (
                "https://raw.githubusercontent.com/innate2adaptive/Decombinator-Tags-FASTAs/master/"
                + expected_file
            )
            urllib.request.urlopen(fl)  # Request URL, see whether is found
            fl = urllib.request.urlretrieve(fl)[0]
        except:
            print(
                "Cannot find following file locally or online:", expected_file
            )
            print(
                "Please either run Decombinator with internet access, or point Decombinator to local copies of the tag and FASTA files with the '-tf' flag."
            )
            sys.exit()
        # Return opened file, for either FASTA or tag file parsing
    return fl


def readfq(fp):
    """
    readfq(file):Heng Li's Python implementation of his readfq function
    https://github.com/lh3/readfq/blob/master/readfq.py
    """
    # i = 1
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in ">@":  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in "@+>":
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != "+":  # this is a fasta record
            yield name, "".join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = "".join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, "".join(seqs)
                    # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


#####################################
############# DECOMBINE #############
#####################################


def vanalysis(read):

    hold_v = v_key.findall(read)

    if hold_v:
        if len(hold_v) > 1:
            counts["multiple_v_matches"] += 1
            return

        v_match = v_seqs.index(hold_v[0][0])  # Assigns VJ
        temp_end_v = (
            hold_v[0][1] + jump_to_end_v[v_match] - 1
        )  # Finds where the end of a full V would be

        v_seq_start = hold_v[0][1]
        end_v_v_dels = get_v_deletions(read, v_match, temp_end_v, v_regions)
        if end_v_v_dels:  # If the number of deletions has been found
            return v_match, end_v_v_dels[0], end_v_v_dels[1], v_seq_start

    else:

        hold_v1 = half1_v_key.findall(read)

        if hold_v1:
            for i in range(len(hold_v1)):
                indices = [
                    y for y, x in enumerate(half1_v_seqs) if x == hold_v1[i][0]
                ]
                for k in indices:
                    if len(v_seqs[k]) == len(
                        read[
                            hold_v1[i][1] : hold_v1[i][1]
                            + len(v_seqs[half1_v_seqs.index(hold_v1[i][0])])
                        ]
                    ):
                        if (
                            lev.hamming(
                                v_seqs[k],
                                read[
                                    hold_v1[i][1] : hold_v1[i][1]
                                    + len(v_seqs[k])
                                ],
                            )
                            <= 1
                        ):
                            counts["verr2"] += 1
                            v_match = k
                            temp_end_v = (
                                hold_v1[i][1] + jump_to_end_v[v_match] - 1
                            )  # Finds where the end of a full V would be
                            end_v_v_dels = get_v_deletions(
                                read, v_match, temp_end_v, v_regions
                            )
                            if end_v_v_dels:
                                v_seq_start = hold_v1[i][1]
                                return (
                                    v_match,
                                    end_v_v_dels[0],
                                    end_v_v_dels[1],
                                    v_seq_start,
                                )
            counts["foundv1notv2"] += 1
            return

        else:

            hold_v2 = half2_v_key.findall(read)
            if hold_v2:
                for i in range(len(hold_v2)):
                    indices = [
                        y
                        for y, x in enumerate(half2_v_seqs)
                        if x == hold_v2[i][0]
                    ]
                    for k in indices:
                        if len(v_seqs[k]) == len(
                            read[
                                hold_v2[i][1]
                                - v_half_split : hold_v2[i][1]
                                - v_half_split
                                + len(
                                    v_seqs[half2_v_seqs.index(hold_v2[i][0])]
                                )
                            ]
                        ):
                            if (
                                lev.hamming(
                                    v_seqs[k],
                                    read[
                                        hold_v2[i][1]
                                        - v_half_split : hold_v2[i][1]
                                        + len(v_seqs[k])
                                        - v_half_split
                                    ],
                                )
                                <= 1
                            ):
                                counts["verr1"] += 1
                                v_match = k
                                temp_end_v = (
                                    hold_v2[i][1]
                                    + jump_to_end_v[v_match]
                                    - v_half_split
                                    - 1
                                )  # Finds where the end of a full V would be
                                end_v_v_dels = get_v_deletions(
                                    read, v_match, temp_end_v, v_regions
                                )
                                if end_v_v_dels:
                                    v_seq_start = hold_v2[i][1] - v_half_split
                                    return (
                                        v_match,
                                        end_v_v_dels[0],
                                        end_v_v_dels[1],
                                        v_seq_start,
                                    )
                counts["foundv2notv1"] += 1
                return

            else:
                counts["no_vtags_found"] += 1
                return


def janalysis(read, end_of_v):

    hold_j = j_key.findall(read)

    if hold_j:
        if len(hold_j) > 1:
            counts["multiple_j_matches"] += 1
            return

        j_match = j_seqs.index(hold_j[0][0])  # Assigns J
        temp_start_j = (
            hold_j[0][1] - jump_to_start_j[j_match]
        )  # Finds where the start of a full J would be

        j_seq_end = hold_j[0][1] + len(hold_j[0][0])

        start_j_j_dels = get_j_deletions(
            read, j_match, temp_start_j, j_regions, end_of_v
        )

        if start_j_j_dels:  # If the number of deletions has been found
            return j_match, start_j_j_dels[0], start_j_j_dels[1], j_seq_end

    else:

        hold_j1 = half1_j_key.findall(read)
        if hold_j1:
            for i in range(len(hold_j1)):
                indices = [
                    y for y, x in enumerate(half1_j_seqs) if x == hold_j1[i][0]
                ]
                for k in indices:
                    if len(j_seqs[k]) == len(
                        read[
                            hold_j1[i][1] : hold_j1[i][1]
                            + len(j_seqs[half1_j_seqs.index(hold_j1[i][0])])
                        ]
                    ):
                        if (
                            lev.hamming(
                                j_seqs[k],
                                read[
                                    hold_j1[i][1] : hold_j1[i][1]
                                    + len(j_seqs[k])
                                ],
                            )
                            <= 1
                        ):
                            counts["jerr2"] += 1
                            j_match = k
                            temp_start_j = (
                                hold_j1[i][1] - jump_to_start_j[j_match]
                            )  # Finds where the start of a full J would be
                            j_seq_end = (
                                hold_j1[i][1]
                                + len(hold_j1[i][0])
                                + j_half_split
                            )
                            start_j_j_dels = get_j_deletions(
                                read,
                                j_match,
                                temp_start_j,
                                j_regions,
                                end_of_v,
                            )
                            if start_j_j_dels:
                                return (
                                    j_match,
                                    start_j_j_dels[0],
                                    start_j_j_dels[1],
                                    j_seq_end,
                                )
            counts["foundj1notj2"] += 1
            return

        else:
            hold_j2 = half2_j_key.findall(read)
            if hold_j2:
                for i in range(len(hold_j2)):
                    indices = [
                        y
                        for y, x in enumerate(half2_j_seqs)
                        if x == hold_j2[i][0]
                    ]
                    for k in indices:
                        if len(j_seqs[k]) == len(
                            read[
                                hold_j2[i][1]
                                - j_half_split : hold_j2[i][1]
                                - j_half_split
                                + len(
                                    j_seqs[half2_j_seqs.index(hold_j2[i][0])]
                                )
                            ]
                        ):
                            if (
                                lev.hamming(
                                    j_seqs[k],
                                    read[
                                        hold_j2[i][1]
                                        - j_half_split : hold_j2[i][1]
                                        + len(j_seqs[k])
                                        - j_half_split
                                    ],
                                )
                                <= 1
                            ):
                                counts["jerr1"] += 1
                                j_match = k
                                temp_start_j = (
                                    hold_j2[i][1]
                                    - jump_to_start_j[j_match]
                                    - j_half_split
                                )  # Finds where the start of a full J would be
                                j_seq_end = hold_j2[i][1] + len(hold_j2[i][0])
                                start_j_j_dels = get_j_deletions(
                                    read,
                                    j_match,
                                    temp_start_j,
                                    j_regions,
                                    end_of_v,
                                )
                                if start_j_j_dels:
                                    return (
                                        j_match,
                                        start_j_j_dels[0],
                                        start_j_j_dels[1],
                                        j_seq_end,
                                    )
                counts["foundv2notv1"] += 1
                return

            else:
                counts["no_j_assigned"] += 1
                return


def dcr(read, inputargs):
    """dcr(read): Core function which checks a read (in the given frame) for a rearranged TCR of the specified chain.
    Returns a list giving: V gene index, J gene index, # deletions in V gene, # deletions in J gene,
      insert sequence (between ends of V and J), inter-tag sequence (for collapsing), and its quality scores
    """
    v_seq_start = 0
    j_seq_end = 0

    vdat = vanalysis(read)

    if not vdat:
        return

    end_of_v = vdat[1] + 1
    jdat = janalysis(read, end_of_v)

    if jdat:

        # Filter out rearrangements with indications they probably represent erroneous sequences
        if (
            "N" in read[vdat[3] : jdat[3]] and inputargs["allowNs"] == False
        ):  # Ambiguous base in inter-tag region
            counts["dcrfilter_intertagN"] += 1
        elif (vdat[3] - jdat[3]) >= inputargs[
            "lenthreshold"
        ]:  # Inter-tag length threshold
            counts["dcrfilter_toolong_intertag"] += 1
        elif (
            vdat[2] > (jump_to_end_v[vdat[0]] - len(v_seqs[vdat[0]]))
            or jdat[2] > jump_to_start_j[jdat[0]]
        ):  # Impossible number of deletions
            counts["dcrfilter_imposs_deletion"] += 1
        elif (vdat[3] + len(v_seqs[vdat[0]])) > (
            jdat[3] + len(j_seqs[jdat[0]])
        ):  # Overlapping tags
            counts["dcrfilter_tag_overlap"] += 1

        else:
            vj_details = [
                vdat[0],
                jdat[0],
                vdat[2],
                jdat[2],
                read[vdat[1] + 1 : jdat[1]],
                vdat[3],
                jdat[3],
            ]
            return vj_details

    else:
        counts["VJ_assignment_failed"] += 1
        return


###########################################################
############# ANCILLARY DECOMBINING FUNCTIONS #############
###########################################################


def import_tcr_info(inputargs):
    """import_tcr_info: Gathers the required TCR chain information for Decombining"""

    # Get chain information
    global chainnams, chain, counts
    counts = coll.Counter()
    chainnams = {"a": "alpha", "b": "beta", "g": "gamma", "d": "delta"}

    # Detect whether chain specified in filename
    inner_filename_chains = [
        x for x in chainnams.values() if x in inputargs["infile"].lower()
    ]
    if len(inner_filename_chains) == 1:
        counts["chain_detected"] = 1

    if inputargs["chain"]:
        if inputargs["chain"].upper() in ["A", "ALPHA", "TRA", "TCRA"]:
            chain = "a"
        elif inputargs["chain"].upper() in ["B", "BETA", "TRB", "TCRB"]:
            chain = "b"
        elif inputargs["chain"].upper() in ["G", "GAMMA", "TRG", "TCRG"]:
            chain = "g"
        elif inputargs["chain"].upper() in ["D", "DELTA", "TRD", "TCRD"]:
            chain = "d"
        else:
            print(nochain_error)
            sys.exit()
    else:

        # If no chain provided, try and infer from filename
        if counts["chain_detected"] == 1:
            chain = inner_filename_chains[0][0]

        else:
            nochain_error = "TCR chain not recognised. \n \
      Please either include (one) chain name in the file name (i.e. alpha/beta/gamma/delta),\n \
      or use the '-c' flag with an explicit chain option (a/b/g/d, case-insensitive)."
            print(nochain_error)
            sys.exit()

    #################################################
    ############# GET GENES, BUILD TRIE #############
    #################################################

    print("Importing TCR", chainnams[chain], "gene sequences...")

    # First check that valid tag/species combinations have been used
    if inputargs["tags"] == "extended" and inputargs["species"] == "mouse":
        print(
            "Please note that there is currently no extended tag set for mouse TCR genes.\n \
    Decombinator will now switch the tag set in use from 'extended' to 'original'.\n \
    In future, consider editing the script to change the default, or use the appropriate flags (-sp mouse -tg original)."
        )
        inputargs["tags"] = "original"

    if inputargs["tags"] == "extended" and (chain == "g" or chain == "d"):
        print(
            "Please note that there is currently no extended tag set for gamma/delta TCR genes.\n \
    Decombinator will now switch the tag set in use from 'extended' to 'original'.\n \
    In future, consider editing the script to change the default, or use the appropriate flags."
        )
        inputargs["tags"] = "original"

    # Set tag split position, and check tag set. Note that original tags use shorter length J half tags, as these tags were originally shorter.
    global v_half_split, j_half_split
    if inputargs["tags"] == "extended":
        v_half_split, j_half_split = [10, 10]
    elif inputargs["tags"] == "original":
        v_half_split, j_half_split = [10, 6]
    else:
        print(
            "Tag set unrecognised; should be either 'extended' or 'original' for human, or just 'original' for mouse. \n \
    Please check tag set and species flag."
        )
        sys.exit()

    # Check species information
    if inputargs["species"] not in ["human", "mouse"]:
        print(
            "Species not recognised. Please select either 'human' (default) or 'mouse'.\n \
    If mouse is required by default, consider changing the default value in the script."
        )
        sys.exit()

    # Look for tag and V/J fasta and tag files: if these cannot be found in the working directory, source them from GitHub repositories
    # Note that fasta/tag files fit the pattern "species_tagset_gene.[fasta/tags]"
    # I.e. "[human/mouse]_[extended/original]_TR[A/B/G/D][V/J].[fasta/tags]"

    for gene in ["v", "j"]:
        # Get FASTA data
        fasta_file = read_tcr_file(
            inputargs["species"],
            inputargs["tags"],
            gene,
            "fasta",
            inputargs["tagfastadir"],
        )
        globals()[gene + "_genes"] = list(SeqIO.parse(fasta_file, "fasta"))

        globals()[gene + "_regions"] = []
        for g in range(0, len(globals()[gene + "_genes"])):
            globals()[gene + "_regions"].append(
                globals()[gene + "_genes"][g].seq.upper()
            )

        # Get tag data
        tag_file = read_tcr_file(
            inputargs["species"],
            inputargs["tags"],
            gene,
            "tags",
            inputargs["tagfastadir"],
        )  # get tag data
        tag_data = open(tag_file, "r")
        if gene == "v":
            jumpfunction = "jump_to_end_v"
        elif gene == "j":
            jumpfunction = "jump_to_start_j"
        (
            globals()[gene + "_seqs"],
            globals()["half1_" + gene + "_seqs"],
            globals()["half2_" + gene + "_seqs"],
            globals()[jumpfunction],
        ) = globals()["get_" + gene + "_tags"](
            tag_data, globals()[gene + "_half_split"]
        )
        tag_data.close()

        # Build Aho-Corasick tries
        globals()[gene + "_builder"] = AcoraBuilder()
        for i in range(0, len(globals()[gene + "_seqs"])):
            globals()[gene + "_builder"].add(
                str(globals()[gene + "_seqs"][i])
            )  # Add all V tags to keyword trie
        globals()[gene + "_key"] = globals()[gene + "_builder"].build()

        # And tries for split, half-tags
        globals()[gene + "_half1_builder"] = AcoraBuilder()
        for i in range(0, len(globals()["half1_" + gene + "_seqs"])):
            globals()[gene + "_half1_builder"].add(
                str(globals()["half1_" + gene + "_seqs"][i])
            )
        globals()["half1_" + gene + "_key"] = globals()[
            gene + "_half1_builder"
        ].build()

        globals()[gene + "_half2_builder"] = AcoraBuilder()
        for i in range(0, len(globals()["half2_" + gene + "_seqs"])):
            globals()[gene + "_half2_builder"].add(
                str(globals()["half2_" + gene + "_seqs"][i])
            )
        globals()["half2_" + gene + "_key"] = globals()[
            gene + "_half2_builder"
        ].build()


def get_v_deletions(read, v_match, temp_end_v, v_regions_cut):
    # This function determines the number of V deletions in sequence read
    # by comparing it to v_match, beginning by making comparisons at the
    # end of v_match and at position temp_end_v in read.
    function_temp_end_v = temp_end_v
    pos = (
        len(v_regions_cut[v_match]) - 10
    )  # changed from -1 for new checking technique
    is_v_match = 0

    # Catch situations in which the temporary end of the V exists beyond the end of the read
    if function_temp_end_v >= len(read):
        counts["v_del_failed_tag_at_end"] += 1
        return

    function_temp_end_v += 1
    num_del = 0

    while is_v_match == 0 and 0 <= function_temp_end_v < len(read):
        # Require a 10 base match to determine where end of germ-line sequence lies
        if (
            str(v_regions_cut[v_match])[pos : pos + 10]
            == read[function_temp_end_v - 10 : function_temp_end_v]
        ):
            is_v_match = 1
            deletions_v = num_del
            end_v = temp_end_v - num_del
        else:
            pos -= 1
            num_del += 1
            function_temp_end_v -= 1

    if is_v_match == 1:
        return [end_v, deletions_v]
    else:
        counts["v_del_failed"] += 1
        return


def get_j_deletions(read, j_match, temp_start_j, j_regions_cut, end_of_v):
    # This function determines the number of J deletions in sequence read
    # by comparing it to j_match, beginning by making comparisons at the
    # end of j_match and at position temp_end_j in read.
    function_temp_start_j = temp_start_j
    pos = 0
    is_j_match = 0
    while is_j_match == 0 and 0 <= function_temp_start_j + 2 < len(str(read)):
        # in the case of no detectable insertions, where nucleotide junctions could be derived from either gene,
        # nucleotides will be deemed to have derived from the V gene, and count towards deletions from J.
        if function_temp_start_j < end_of_v:
            pos += 1
            function_temp_start_j += 1
        # Require a 10 base match to determine where end of germ-line sequence lies
        elif (
            str(j_regions_cut[j_match])[pos : pos + 10]
            == read[function_temp_start_j : function_temp_start_j + 10]
        ):
            is_j_match = 1
            deletions_j = pos
            start_j = function_temp_start_j
        else:
            pos += 1
            function_temp_start_j += 1

    if is_j_match == 1:
        return [start_j, deletions_j]
    else:
        counts["j_del_failed"] += 1
        return


def get_v_tags(file_v, half_split):
    # """Read V tags in from file"""
    v_seqs = []  # Holds all V tags
    jump_to_end_v = (
        []
    )  # Holds the number of jumps to make to look for deletions for each V region once the corresponding tag has been found
    for line in file_v:
        elements = line.rstrip(
            "\n"
        )  # Turns every element in a text file line separated by a space into elements in a list
        v_seqs.append(
            elements.split()[0]
        )  # Adds elements in first column iteratively
        jump_to_end_v.append(
            int(elements.split()[1])
        )  # Adds elements in second column iteratively

    half1_v_seqs = []
    half2_v_seqs = []

    for i in range(len(v_seqs)):
        half1_v_seqs.append(v_seqs[i][0:half_split])
        half2_v_seqs.append(v_seqs[i][half_split:])

    return [v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v]


def get_j_tags(file_j, half_split):
    """Read J tags in from file"""
    j_seqs = []  # Holds all J tags
    jump_to_start_j = (
        []
    )  # Holds the number of jumps to make to look for deletions for each J region once the corresponding tag has been found

    for line in file_j:
        elements = line.rstrip("\n")
        j_seqs.append(elements.split()[0])
        jump_to_start_j.append(int(elements.split()[1]))

    half1_j_seqs = []
    half2_j_seqs = []

    for j in range(len(j_seqs)):
        half1_j_seqs.append(j_seqs[j][0:half_split])
        half2_j_seqs.append(j_seqs[j][half_split:])

    return [j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j]


def sort_permissions(fl):
    # Need to ensure proper file permissions on output data
    # If users are running pipeline through Docker might otherwise require root access
    if oct(os.stat(fl).st_mode)[4:] != "666":
        os.chmod(fl, 0o666)


##########################################################
############# READ IN COMMAND LINE ARGUMENTS #############
##########################################################


def decombinator(inputargs: dict) -> list:
    """Function wrapper for decombinator."""

    print("Running Decombinator version", metadata.version("decombinator"))

    opener = opener_check(inputargs)

    # Get TCR gene information
    import_tcr_info(inputargs)
    samplenam = str(inputargs["infile"].split(".")[0])

    if (
        os.sep in samplenam
    ):  # Cope with situation where specified FQ file is in a subdirectory
        samplenam = samplenam.split(os.sep)[-1]

    if inputargs["suppresssummary"] == False:

        logpath = inputargs["outpath"] + f"Logs{os.sep}"

        # Check for directory and make summary file
        if not os.path.exists(logpath):
            os.makedirs(logpath)
        date = strftime("%Y_%m_%d")

        # Check for existing date-stamped file
        summaryname = logpath + date + "_"
        if inputargs["chain"]:
            summaryname += chainnams[chain] + "_"
        summaryname += samplenam + "_Decombinator_Summary.csv"

    # Brief FASTQ sanity check
    if inputargs["dontcheck"] == False:
        fastq_check(inputargs, opener, samplenam, summaryname, logpath)

    # Get Barcode length
    bclength = inputargs["bclength"]

    counts["start_time"] = time()

    #########################################################
    ############# SCROLL THROUGH FILE & ANALYSE #############
    #########################################################

    print("Decombining FASTQ data...")

    suffix = "." + inputargs["extension"]

    # If chain had not been autodetected, write it out into output file
    if counts["chain_detected"] == 1:
        name_results = inputargs["prefix"] + samplenam
    else:
        name_results = inputargs["prefix"] + samplenam + "_" + chainnams[chain]

    if inputargs["nobarcoding"] == False:
        stemplate = string.Template(
            "$v $j $del_v $del_j $nt_insert $seqid $tcr_seq $tcr_qual $barcode $barqual"
        )
    else:
        stemplate = string.Template("$v $j $del_v $del_j $nt_insert")
    found_tcrs = coll.Counter()
    ############################################################################################
    ##################################################################################################

    ##################################################################################
    # Scroll through input file and find TCRs

    outdata = []
    start_time = time()
    if inputargs["nobarcoding"] == False:
        if inputargs["bc_read"] == "R2":
            fq1 = readfq(opener(inputargs["infile"], "rt"))
            fq2 = readfq(
                opener(inputargs["infile"].replace("1.f", "2.f"), "rt")
            )
        elif inputargs["bc_read"] == "R1":
            fq1 = readfq(opener(inputargs["infile"], "rt"))
            fq2 = fq1

        fqs = (fq1, fq2)
        zipfqs = zip(fq1, fq2)

        for records in zipfqs:
            record1, record2 = records
            if inputargs["bc_read"] == "R2":
                readid = record1[0]
                vdj = record1[1]
                vdjqual = record1[2]
                bc = record2[1][:bclength]
                bcQ = record2[2][:bclength]

            elif inputargs["bc_read"] == "R1":
                readid = record1[0]
                vdj = record1[1][bclength:]
                vdjqual = record1[2][bclength:]
                bc = record1[1][0:bclength]
                bcQ = record1[2][0:bclength]

            if inputargs["nobarcoding"] == False:
                if (
                    "N" in bc and inputargs["allowNs"] == False
                ):  # Ambiguous base in barcode region
                    counts["dcrfilter_barcodeN"] += 1

            counts["read_count"] += 1
            if (
                counts["read_count"] % 100000 == 0
                and inputargs["dontcount"] == False
            ):
                print("\t read", counts["read_count"])

            # Get details of the VJ recombination
            if inputargs["orientation"] == "reverse":
                recom = dcr(revcomp(vdj), inputargs)
                frame = "reverse"
            elif inputargs["orientation"] == "forward":
                recom = dcr(vdj, inputargs)
                frame = "forward"
            elif inputargs["orientation"] == "both":
                recom = dcr(revcomp(vdj), inputargs)
                frame = "reverse"
                if not recom:
                    recom = dcr(vdj, inputargs)
                    frame = "forward"

            if recom:
                counts["vj_count"] += 1

                if frame == "reverse":
                    tcrseq = revcomp(vdj)[recom[5] : recom[6]]
                    tcrQ = vdjqual[::-1][recom[5] : recom[6]]
                elif frame == "forward":
                    tcrseq = vdj[recom[5] : recom[6]]
                    tcrQ = vdjqual[recom[5] : recom[6]]

                if inputargs["nobarcoding"] == False:
                    dcr_output = [
                        str(recom[0]),
                        str(recom[1]),
                        str(recom[2]),
                        str(recom[3]),
                        recom[4],
                        readid,
                        tcrseq,
                        tcrQ,
                        bc,
                        bcQ,
                    ]
                    outdata.append(dcr_output)

                else:  # TODO: create non-barcode alternative OR delete non-barcode methodology
                    dcr_string = stemplate.substitute(
                        v=str(recom[0]) + ",",
                        j=str(recom[1]) + ",",
                        del_v=str(recom[2]) + ",",
                        del_j=str(recom[3]) + ",",
                        nt_insert=recom[4],
                    )
                    found_tcrs[dcr_string] += 1
                    outdata.append(dcr_string)

    if inputargs["nobarcoding"] == True:
        # Write out non-barcoded results, with frequencies
        if inputargs["extension"] == "n12":
            print(
                "Non-barcoding option selected, but default output file extension (n12) detected. Automatically changing to 'nbc'."
            )
            suffix = ".nbc"
        for x in found_tcrs.most_common():
            outdata.append(x[0] + ", " + str(found_tcrs[x[0]]))

    counts["end_time"] = time()
    timetaken = counts["end_time"] - counts["start_time"]

    ##############################################
    ############# WRITE SUMMARY DATA #############
    ##############################################

    print(
        "Analysed",
        "{:,}".format(counts["read_count"]),
        "reads, finding",
        "{:,}".format(counts["vj_count"]),
        chainnams[chain],
        "VJ rearrangements",
    )
    print("Reading from", inputargs["infile"] + ", writing to variable")
    print("Took", str(round(timetaken, 2)), "seconds")

    # Write data to summary file
    if inputargs["suppresssummary"] == False:
        if not os.path.exists(summaryname):
            summaryfile = open(summaryname, "wt")
        else:
            # If one exists, start an incremental day stamp
            for i in range(2, 10000):
                summaryname = logpath + date + "_"
                if inputargs["chain"]:
                    summaryname += chainnams[chain] + "_"
                summaryname += (
                    samplenam + "_Decombinator_Summary" + str(i) + ".csv"
                )
                if not os.path.exists(summaryname):
                    summaryfile = open(summaryname, "wt")
                    break

        inout_name = (
            "_".join(f"{samplenam}".split("_")[:-1]) + f"_{chainnams[chain]}"
        )

        # Generate string to write to summary file
        summstr = (
            "Property,Value\nDirectory,"
            + os.getcwd()
            + "\nInputFile,"
            + inout_name
            + "\nOutputFile,"
            + inout_name
            + "\nDateFinished,"
            + date
            + "\nTimeFinished,"
            + strftime("%H:%M:%S")
            + "\nTimeTaken(Seconds),"
            + str(round(timetaken, 2))
            + "\n\nInputArguments:,\n"
        )
        for s in [
            "species",
            "chain",
            "extension",
            "tags",
            "dontgzip",
            "allowNs",
            "orientation",
            "lenthreshold",
            "bc_read",
            "bclength",
        ]:
            summstr = summstr + s + "," + str(inputargs[s]) + "\n"

        counts["pc_decombined"] = counts["vj_count"] / counts["read_count"]

        summstr = (
            summstr
            + "\nNumberReadsInput,"
            + str(counts["read_count"])
            + "\nNumberReadsDecombined,"
            + str(counts["vj_count"])
            + "\nPercentReadsDecombined,"
            + str(round(counts["pc_decombined"], 3))
        )

        # Half tag matching details
        summstr = (
            summstr
            + "\n\nReadsAssignedUsingHalfTags:,\nV1error,"
            + str(counts["verr1"])
            + "\nV2error,"
            + str(counts["verr2"])
            + "\nJ1error,"
            + str(counts["jerr1"])
            + "\nJ2error,"
            + str(counts["jerr2"])
        )

        # Number reads filtered out
        summstr = (
            summstr
            + "\n\nReadsFilteredOut:,\nAmbiguousBaseCall(DCR),"
            + str(counts["dcrfilter_intertagN"])
            + "\nAmbiguousBaseCall(Barcode),"
            + str(counts["dcrfilter_barcodeN"])
            + "\nOverlongInterTagSeq,"
            + str(counts["dcrfilter_toolong_intertag"])
            + "\nImpossibleDeletions,"
            + str(counts["dcrfilter_imposs_deletion"])
            + "\nOverlappingTagBoundaries,"
            + str(counts["dcrfilter_tag_overlap"])
        )
        ##########################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##########################
        summstr = (
            summstr
            + "\n\nReadsFailedAssignment:,\nMultipleVtagMatches,"
            + str(counts["multiple_v_matches"])
            + "\nVTagAtEndRead,"
            + str(counts["v_del_failed_tag_at_end"])
            + "\nVDeletionsUndetermined,"
            + str(counts["v_del_failed"])
            + "\nFoundV1HalfTagNotV2,"
            + str(counts["foundv1notv2"])
            + "\nFoundV2HalfTagNotV1,"
            + str(counts["foundv2notv1"])
            + "\nNoVDetected,"
            + str(counts["no_vtags_found"])
            + "\nMultipleJTagMatches,"
            + str(counts["multiple_j_matches"])
            + "\nJDeletionsUndermined,"
            + str(counts["j_del_failed"])
            + "\nFoundJ1HalfTagNotJ2,"
            + str(counts["foundj1notj2"])
            + "\nFoundJ2HalfTagNotJ1,"
            + str(counts["foundj2notj1"])
            + "\nNoJDetected,"
            + str(counts["no_j_assigned"])
        )
        # + "\nVJGeneAssignmentFailed," + str(counts['VJ_assignment_failed'])

        print(summstr, file=summaryfile)
        summaryfile.close()
        sort_permissions(summaryname)

    return outdata


if __name__ == "__main__":
    print(
        "Calling Decombinator from the shell has been depreciated as of Decombinator V4.3. \
        Please check the README on how to update your script, or alternativley change branch to decombinator_v4.2 \
        which retains this functionality."
    )
