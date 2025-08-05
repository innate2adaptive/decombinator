![Decombinator Logo](./assets/decombinator-logo.png)

# `decombinator`

[![License](https://img.shields.io/badge/license-MIT-blue)](https://raw.githubusercontent.com/innate2adaptive/Decombinator/master/LICENCE)
<!-- [![Latest release](https://img.shields.io/pypi/v/decombinator)](https://pypi.org/p/decombinator) -->
<!-- [![Python versions](https://img.shields.io/pypi/pyversions/decombinator)](https://pypi.org/project/decombinator/) -->
<!-- ![Tests](https://github.com/innate2adaptive/decombinator/.github/workflows/package-build-test.yml/badge.svg) -->

**Innate2Adaptive lab @ University College London, 2024**

**Written by (in alphabetical order) Benny Chain, James M. Heather, Katharine Best, Matthew V. Cowley, Mazlina Ismail, Niclas Thomas, Tahel Ronel, Theres Oakes, and Thomas Peacock.**

`decombinator` is a fast and efficient tool for the analysis of T-cell receptor (TCR) repertoire sequences produced by deep sequencing.
It accepts TCR repertoire sequencing (TCRseq) FASTQ data and returns [AIRR](https://docs.airr-community.org/en/stable/datarep/rearrangements.html) compliant files detailing TCR V(D)J recombination and counts.

TCRseq offers a powerful means to investigate biological samples to observe sequence distributions.
However current high-throughput sequencing (HTS) technologies can produce large amounts of raw data which will also unavoidably contain errors relative to the original input molecules.
`decombinator` addresses the problem of large datasets through speed - employing a rapid and [highly efficient string matching algorithm](https://figshare.com/articles/Aho_Corasick_String_Matching_Video/771968) to search the FASTQ files produced by HTS machines for rearranged TCR sequences.
The central algorithm searches for 'tag' sequences, the presence of which uniquely indicates the inclusion of particular V or J genes in a recombination.
If V and J tags are found, `decombinator` can then deduce where the ends of the germline V and J gene sections are (i.e. how much nucleotide removal occurred during V(D)J recombination), and what nucleotide sequence (the 'insert sequence') remains between the two.
These five pieces of information - the V and J genes used, how many deletions each had and the insert sequence - contain all of the information required to reconstruct the whole TCR nucleotide sequence, in a more readily stored and analysed way.
The pipeline also handles FASTQ data with many reads of the same molecule due to over-amplification.
It does this by using unique molecular identifiers to "collapse" down duplicate reads.

Running `decombinator` is easy. All that is required is HTS read data and a few arguments specifying the chemistry of your prepared library:

```shell
decombinator pipeline -in XXXX.fq -c b -br R2 -bl 42 -ol M13
```

## Installation

To install `decombinator` and all required packages, simply create a fresh virtual environment, activate it, and run:

```shell
pip install decombinator
```

`decombinator` also requires a number of additional files, which contain information regarding the V and J gene sequences, their tags, and the locations and motifs which define their CDR3 regions. By default `decombinator` downloads these files from [the git repository where they are maintained](https://github.com/innate2adaptive/Decombinator-Tags-FASTAs), which requires a working internet connection. In order to run `decombinator` offline, these files must be downloaded to a local location, and either stored within the directory where you wish to run `decombinator` or specified using the appropriate command line flag (`-tfdir`). These additional files can also be downloaded via `git clone`:

```bash
git clone https://github.com/innate2adaptive/Decombinator-Tags-FASTAs.git
```

The current version of `decombinator` has tag sets for the analysis of alpha/beta and gamma/delta TCR repertoires from both mouse and man. In addition to the original tag set developed for the 2013 Thomas *et al* Bioinformatics paper, an 'extended' tag set has been developed for human alpha/beta repertoires, which covers 'non-functional' V and J genes (ORFs and pseudogenes, according to IMGT nomenclature) as well as just the functional.

## Basic Usage

### `pipeline`

`decombinator` can be run in 4 different modes depending on your use case. Firstly, `pipeline` mode can be used via:

```shell
decombinator pipeline -in XXXX.fq -c b -br R2 -bl 42 -ol M13
```

<details>
  <summary><code>decombinator pipeline</code> CLI arguments</summary>
  
  | Option                               | Description                                                                                                                                                  |
|--------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-h`, `--help`                       | Show this help message and exit                                                                                                                               |
| `-s`, `--suppresssummary`            | Suppress the production of summary data log/file                                                                                                              |
| `-dz`, `--dontgzip`                  | Stop the output FASTQ files automatically being compressed with gzip                                                                                          |
| `-dc`, `--dontcount`                 | Stop/Block printing the running count                                                                                                                         |
| `-op OUTPATH`, `--outpath OUTPATH`   | Path to output directory, writes to directory script was called in by default                                                                                 |
| `-c CHAIN`, `--chain CHAIN`          | TCR chain (a/b/g/d)                                                                                                                                           |
| `-pf PREFIX`, `--prefix PREFIX`      | Specify the prefix of the output DCR file. Default = "dcr_"                                                                                                    |
| `-in INFILE`, `--infile INFILE`      | Correctly demultiplexed/processed FASTQ file containing TCR reads                                                                                              |
| `-br BC_READ`, `--bc_read BC_READ`   | Which read has bar code (R1,R2). If used, ensure read selected is present in the same directory as the file specified by -in.                                                                                                                                |
| `-dk`, `--dontcheck`                 | Skip the FASTQ check                                                                                                                                           |
| `-ex EXTENSION`, `--extension EXTENSION` | Specify the file extension of the output DCR file. Default = "n12"                                                                                         |
| `-or ORIENTATION`, `--orientation ORIENTATION` | Specify the orientation to search in (forward/reverse/both). Default = reverse                                                                  |
| `-tg TAGS`, `--tags TAGS`            | Specify which Decombinator tag set to use (extended or original). Default = extended                                                                           |
| `-sp SPECIES`, `--species SPECIES`   | Specify which species TCR repertoire the data consists of (human or mouse). Default = human                                                                    |
| `-N`, `--allowNs`                    | Whether to allow VJ rearrangements containing ambiguous base calls ('N'). Default = False                                                                      |
| `-ln LENTHRESHOLD`, `--lenthreshold LENTHRESHOLD` | Acceptable threshold for inter-tag (V to J) sequence length. Default = 130                                                                      |
| `-tfdir TAGFASTADIR`, `--tagfastadir TAGFASTADIR` | Path to folder containing TCR FASTA and Decombinator tag files, for offline analysis. Default = "Decombinator-Tags-FASTAs".                      |
| `-nbc`, `--nobarcoding`              | Option to run Decombinator without barcoding, i.e. so as to run on data produced by any protocol.                                                              |
| `-bl BCLENGTH`, `--bclength BCLENGTH` | Length of barcode sequence, if applicable. Default is set to 42 bp.                                                                                          |
| `-mq MINBCQ`, `--minbcQ MINBCQ`      | Minimum quality score that barcode nucleotides should be to for that rearrangement to be retained. Default = 20.                                               |
| `-bm BCQBELOWMIN`, `--bcQbelowmin BCQBELOWMIN` | Number of nucleotides per barcode whose quality score are allowed to be below -mq and still be retained. Default = 1.                                 |
| `-aq AVGQTHRESHOLD`, `--avgQthreshold AVGQTHRESHOLD` | Average quality threshold that barcode sequences must remain above for rearrangements to be retained. Default = 30                                  |
| `-lv PERCENTLEVDIST`, `--percentlevdist PERCENTLEVDIST` | Percentage Levenshtein distance that is allowed to estimate whether two sequences within a barcode are derived from the same originator molecule. Default = 10 |
| `-bc BCTHRESHOLD`, `--bcthreshold BCTHRESHOLD` | Number of sequence edits that are allowed to consider two barcodes to be derived from same originator during clustering. Default = 2.                    |
| `-di`, `--dontcheckinput`            | Override the input file sanity check                                                                                                                           |
| `-bd`, `--barcodeduplication`        | Optionally output a file containing the final list of clustered barcodes, and their frequencies                                                                |
| `-pb`, `--positionalbarcodes`        | Instead of inferring random barcode sequences from their context relative to spacer sequences, just take the sequence at the default positions. Useful to salvage runs when R2 quality is terrible. |
| `-ol OLIGO`, `--oligo OLIGO`         | Choose experimental oligo for correct identification of spacers ["M13", "I8", "I8_single", "NEBIO", "TAKARA"] (default: M13)                                                         |
| `-wc`, `--writeclusters`             | Write cluster data to separate cluster files                                                                                                                   |
| `-uh`, `--UMIhistogram`              | Creates histogram of average UMI cluster sizes                                                                                                                 |
| `-npf`, `--nonproductivefilter`      | Filter out non-productive reads from the output                                                                                                                |


</details>

In this mode, all three main components of the pipeline are applied to the data sequentially: `decombine`, `collapse`, and `translate`.
For the full details of each of these components please see the respective section below.

If instead, you wish to run just one component of the pipeline, these are all directly accessible via their respective sub-parser.

### `decombine`

```shell
decombinator decombine -in XXXX.fq -c b -br R2 -bl 42
```

<details>
  <summary><code>decombinator decombine</code> CLI arguments</summary>

| Option                               | Description                                                                                                                                                  |
|--------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-h`, `--help`                       | Show this help message and exit                                                                                                                               |
| `-s`, `--suppresssummary`            | Suppress the production of summary data log/file                                                                                                              |
| `-dz`, `--dontgzip`                  | Stop the output FASTQ files automatically being compressed with gzip                                                                                          |
| `-dc`, `--dontcount`                 | Stop/Block printing the running count                                                                                                                         |
| `-op OUTPATH`, `--outpath OUTPATH`   | Path to output directory, writes to directory script was called in by default                                                                                 |
| `-c CHAIN`, `--chain CHAIN`          | TCR chain (a/b/g/d)                                                                                                                                           |
| `-pf PREFIX`, `--prefix PREFIX`      | Specify the prefix of the output DCR file. Default = "dcr_"                                                                                                    |
| `-in INFILE`, `--infile INFILE`      | Correctly demultiplexed/processed FASTQ file containing TCR reads                                                                                              |
| `-br BC_READ`, `--bc_read BC_READ`   | Which read has bar code (R1,R2). If used, ensure read selected is present in the same directory as the file specified by -in.                                                                                                                                |
| `-dk`, `--dontcheck`                 | Skip the FASTQ check                                                                                                                                           |
| `-ex EXTENSION`, `--extension EXTENSION` | Specify the file extension of the output DCR file. Default = "n12"                                                                                         |
| `-or ORIENTATION`, `--orientation ORIENTATION` | Specify the orientation to search in (forward/reverse/both). Default = reverse                                                                  |
| `-tg TAGS`, `--tags TAGS`            | Specify which Decombinator tag set to use (extended or original). Default = extended                                                                           |
| `-sp SPECIES`, `--species SPECIES`   | Specify which species TCR repertoire the data consists of (human or mouse). Default = human                                                                    |
| `-N`, `--allowNs`                    | Whether to allow VJ rearrangements containing ambiguous base calls ('N'). Default = False                                                                      |
| `-ln LENTHRESHOLD`, `--lenthreshold LENTHRESHOLD` | Acceptable threshold for inter-tag (V to J) sequence length. Default = 130                                                                      |
| `-tfdir TAGFASTADIR`, `--tagfastadir TAGFASTADIR` | Path to folder containing TCR FASTA and Decombinator tag files, for offline analysis. Default = "Decombinator-Tags-FASTAs".                      |
| `-nbc`, `--nobarcoding`              | Option to run Decombinator without barcoding, i.e. so as to run on data produced by any protocol.                                                              |
| `-bl BCLENGTH`, `--bclength BCLENGTH` | Length of barcode sequence, if applicable. Default is set to 42 bp.                                                                                          |

</details>

<details>
  <summary><code>decombine</code> details</summary>

This function performs the key computation of the pipeline, as it searches through demultiplexed reads for rearranged TCR sequences. It looks for short 'tag' sequences (using Aho-Corasick string matching): the presence of a tag uniquely identifies a particular V or J gene. If it finds both a V and a J tag (and the read passes various filters), it assigns the read as recombined, and outputs a five-part Decombinator index (or 'DCR'), which uniquely represents a given TCR rearrangement, plus some addtional UMI related information.

All DCR-containing output files are comma-delimited, with the fields of that five-part classifier containing, in order:
* The V index (which V gene was used)
* The J index
* Number of V deletions (relative to germline)
* Number of J deletions
* Insert sequence (the nucleotide sequence between the end of deleted V and J)

The V and J indices are arbitrary numbers based on the order of the tag sequences in the relevant tag file (using Python indexing, which starts at 0 rather than 1). Also, note that the number of V and J deletions just represents how many bases have been removed from the end of that particular germline gene (as given in the germline FASTA files in the additional file repo); it is entirely possible that more bases were deleted, and just that the same bases have been re-added.
Additionally, there are low frequencies of (predominantly alpha chain) recombinations where there is no detectable insertion, and where the nucleotides at the junction between the germline V and J genes could have derived from either. In such circumstances, the nucleotides will arbitrarily be deemed to have derived from the V gene, and thus count towards deletions from the J, however it is impossible to know which gene originally contributed these residues.

Various additional fields may follow the five-part classifier, but the DCR will always occupy the first five positions. An example identifier, from a human alpha chain file, might look like this:

```bash
1, 22, 9, 0, CTCTA
```

Which corresponds to a rearrangement between TRAV1-2 (V index **1**, with **9** nucleotides deleted) and TRAJ33 (J index **22**, with **0** deletions), with an insert sequence (i.e. non-templated additions to the V and/or the J gene) of '**CTCTA**'. For beta chains, the insert sequence will contain any residual TRBD nucleotides, although as these genes are very short, homologous, and typically highly 'nibbled', they are often impossible to differentiate.

Produces a list of lists with the following entries:
  1. V index
  2. J index
  3. Number of V deletions
  4. Number of J deletions
  5. insert
  6. ID
  7. TCR sequence
  8. TCR quality
  9. barcode sequence
  10. barcode quality

**  NB The TCR sequence given here is the 'inter-tag' region, i.e. the sequence between the start of the found V tag and the end of the found J tag.

</details>

### `collapse`

```shell
decombinator collapse -in XXXX.n12 -c b -ol M13
```

<details>
  <summary><code>decombinator collapse</code> CLI arguments</summary>

| Option                               | Description                                                                                                                                                  |
|--------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-h`, `--help`                       | Show this help message and exit                                                                                                                               |
| `-s`, `--suppresssummary`            | Suppress the production of summary data log/file                                                                                                              |
| `-dz`, `--dontgzip`                  | Stop the output FASTQ files automatically being compressed with gzip                                                                                          |
| `-dc`, `--dontcount`                 | Stop/Block printing the running count                                                                                                                         |
| `-op OUTPATH`, `--outpath OUTPATH`   | Path to output directory, writes to directory script was called in by default                                                                                 |
| `-c CHAIN`, `--chain CHAIN`          | TCR chain (a/b/g/d)                                                                                                                                           |
| `-pf PREFIX`, `--prefix PREFIX`      | Specify the prefix of the output DCR file. Default = "dcr_"                                                                                                    |
| `-in INFILE`, `--infile INFILE`      | File containing raw verbose Decombinator output, i.e. 5 part classifier plus barcode and inter-tag sequence and quality strings                                |
| `-mq MINBCQ`, `--minbcQ MINBCQ`      | Minimum quality score that barcode nucleotides should be to for that rearrangement to be retained. Default = 20                                                |
| `-bm BCQBELOWMIN`, `--bcQbelowmin BCQBELOWMIN` | Number of nucleotides per barcode whose quality score are allowed to be below -mq and still be retained. Default = 1                                |
| `-aq AVGQTHRESHOLD`, `--avgQthreshold AVGQTHRESHOLD` | Average quality threshold that barcode sequences must remain above for rearrangements to be retained. Default = 30                                 |
| `-lv PERCENTLEVDIST`, `--percentlevdist PERCENTLEVDIST` | Percentage Levenshtein distance that is allowed to estimate whether two sequences within a barcode are derived from the same originator molecule. Default = 10 |
| `-bc BCTHRESHOLD`, `--bcthreshold BCTHRESHOLD` | Number of sequence edits that are allowed to consider two barcodes to be derived from same originator during clustering. Default = 2                    |
| `-ex EXTENSION`, `--extension EXTENSION` | Specify the file extension of the output DCR file. Default = 'freq'                                                                                         |
| `-N`, `--allowNs`                    | Used to allow VJ rearrangements containing ambiguous base calls ('N')                                                                                          |
| `-ln LENTHRESHOLD`, `--lenthreshold LENTHRESHOLD` | Acceptable threshold for inter-tag (V to J) sequence length                                                                                                  |
| `-di`, `--dontcheckinput`            | Override the input file sanity check                                                                                                                           |
| `-bd`, `--barcodeduplication`        | Optionally output a file containing the final list of clustered barcodes, and their frequencies                                                                |
| `-pb`, `--positionalbarcodes`        | Instead of inferring random barcode sequences from their context relative to spacer sequences, just take the sequence at the default positions. Useful to salvage runs when R2 quality is terrible. |
| `-ol OLIGO`, `--oligo OLIGO`         | Choose experimental oligo for correct identification of spacers ["M13", "I8", "I8_single", "NEBIO", "TAKARA"] (default: M13)                                                         |
| `-wc`, `--writeclusters`             | Write cluster data to separate cluster files                                                                                                                   |
| `-uh`, `--UMIhistogram`              | Creates histogram of average UMI cluster sizes                                                                                                                 |


</details>

<details>
  <summary><code>collapse</code> details</summary>

Takes the output files of Decombinator (run using the barcoding option) and performs collapsing and error correction. This version is a modified version of KB's script collapsinator_20141126.py (That was itself an improved version of the CollapseTCRs.py script used in the Heather et al HIV TCR paper (DOI: 10.3389/fimmu.2015.00644))*
**  Version 4.0.2 includes improved clustering routines measuring the similarity in both barcode and TCR sequence of TCR repertoire data
  
**  NOTE - from version 4.2 this optionally looks for barcode 6NI86N at the beginning of the read; instead of M13_6N_I8_6N_I8
  (i.e. only one spacer).
  This makes it compatible with the multiplex protocol in which the barcode is incorporated in the RT step and is found at the beginning of R1. ** From version V4.2  there is a required additional command line parameter -ol (see below for allowed inputs)
The barcode sequence is contained in one of the additional fields output by `decombine`, which contains the first 42 bases of R2 (if `-br R2 -bl 42` specified). As Illumina sequencing is particularly error-prone in the reverse read, and that reads can be phased (i.e. they do not always begin with the next nucleotide that follows the sequencing primer) our protocol uses known spacer sequences to border the random barcode bases, so that we can identify the actual random bases. The hexameric barcode locations (N6) are determined in reference to the two spacer sequences like so:

```
I8 (spacer) – N6 – I8 – N6 – 2 base overflow (n)
GTCGTGATNNNNNNGTCGTGATNNNNNNnn
```

The collapsing script uses the spacer sequences to identify the exact position of the barcode sequences.

`collapse` performs the following procedures:
* Scrolls through each line of the input object containing DCR, barcode and sequence data.
* Removes TCR reads with forbidden errors, e.g. ambiguous base calls (with user input parameters provided to modify strictness).
* Groups input reads by barcode. Reads with identical barcodes and equivalent inter-tag sequences are grouped together. Equivalence is defined as the Levenshtein distance between two sequences being lower than a given threshold, weighted by the lengths of the compared sequences. Reads with identical barcodes but non-equivalent sequences are grouped separately.
* Each group is assigned the most common inter-tag sequence/DCR combination as the 'true' TCR, as errors are likely to occur during later PCR cycles, and thus will most often be minority variants (see [Bolotin *et al.*, 2012](http://dx.doi.org/10.1002/eji.201242517)).

After this initial grouping, the script estimates the true cDNA frequency. UMIs that are both similar and are associated with a similar TCR are likely to be amplified from the same original DNA molecule and to differ only due to PCR or sequencing error. Consequently, groups with similar barcodes and sequences are then clustered via the following procedure:
* The barcode of each group is compared to the barcode of every other group.
* The expected distribution of distances between UMIs can be modelled as a binomial distribution. Experimentation with simulated datasets found the best threshold for allowing two barcodes to be considered equivalent is when they have a Levenshtein distance of less than 3; a value of 2 is set by default. This can be modified through the user input parameter ```-bc```.
* Groups with barcodes that meet this threshold criteria have their inter-tag sequences compared. Those with equivalent sequences are clustered together. Sequence equivalence is here taken to mean that the two sequences have a Levenshtein distance less than or equal to 10% of the length of the shorter of the two sequences. This percentage can be modified through the user input parameter ```-lv```.
* Upon this merging of groups, the most common inter-tag sequence of the cluster is reassessed and taken as the 'true' TCR. 

Finally, the clusters are collapsed to give the abundance of each TCR in the biological sample.
* A TCR abundance count is calculated for each TCR by counting the number of clusters that have the same sequence but different barcodes (thus representing the same rearrangement originating from multiple input DNA molecules).
* An average UMI count is calculated for each TCR by summing the number of members in each cluster associated with the TCR sequence and dividing by the number of those clusters. This gives a measure that can be used to estimate the robustness of the data for that particular sequence.

Collapsinator outputs 7 fields: the 5-part DCR identifier, the corrected abundance of that TCR in the sample, and the average UMI count for that TCR

</details>

### `translate`

```shell
decombinator translate -in XXXX.freq -c b -ol M13
```

<details>
  <summary><code>decombinator translate</code> CLI arguments</summary>

| Option                               | Description                                                                                                                                                  |
|--------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-h`, `--help`                       | Show this help message and exit                                                                                                                               |
| `-s`, `--suppresssummary`            | Suppress the production of summary data log/file                                                                                                              |
| `-dz`, `--dontgzip`                  | Stop the output FASTQ files automatically being compressed with gzip                                                                                          |
| `-dc`, `--dontcount`                 | Stop/Block printing the running count                                                                                                                         |
| `-op OUTPATH`, `--outpath OUTPATH`   | Path to output directory, writes to directory script was called in by default                                                                                 |
| `-c CHAIN`, `--chain CHAIN`          | TCR chain (a/b/g/d)                                                                                                                                           |
| `-pf PREFIX`, `--prefix PREFIX`      | Specify the prefix of the output DCR file. Default = "dcr_"                                                                                                    |
| `-in INFILE`, `--infile INFILE`      | File containing 5 part classifier plus barcode and inter-tag sequence and quality strings                                                                      |
| `-sp SPECIES`, `--species SPECIES`   | Specify which species TCR repertoire the data consists of (human or mouse). Default = human                                                                    |
| `-tg TAGS`, `--tags TAGS`            | Specify which Decombinator tag set to use (extended or original). Default = extended                                                                           |
| `-npf`, `--nonproductivefilter`      | Filter out non-productive reads from the output                                                                                                                |
| `-tfdir TAGFASTADIR`, `--tagfastadir TAGFASTADIR` | Path to folder containing TCR FASTA and Decombinator tag files, for offline analysis. Default = 'Decombinator-Tags-FASTAs'                           |
| `-nbc`, `--nobarcoding`              | Option to run CD3translator without barcoding, i.e. so as to run on data produced by any protocol.                                                             |


</details>

<details>
  <summary><code>translate</code> details</summary>

This step outputs `.tsv` files in the form of:

| Field | Description | 
|:---:|---|
| sequence_id | A unique identifier for a given rearrangement |
| v_call | V gene used (or multiple, if they cannot be distinguished, comma delimited) |
| d_call | Blank required field (mostly cannot be assigned for TCRb, and rarely useful even then) |
| j_call | J gene used |
| junction_aa | CDR3 junction amino acid sequence |
| duplicate_count | Rearrangement abundance, from Collapsinator |
| sequence | Inferred full-length variable domain nucleotide sequence |
| junction | CDR3 junction nucleotide sequence |
| decombinator_id | Five-field Decombinator identifier |
| rev_comp | Whether rearrangements are reverse complemented (T/F) - this is never the case post-Decombining |
| productive | Whether rearrangement is potentially productive (T/F) |
| sequence_aa | Inferred full-length variable domain amino acid sequence |
| cdr1_aa | Amino acid sequence of CDR1 of the used V gene |
| cdr2_aa | Amino acid sequence of CDR2 of the used V gene |
| vj_in_frame | Whether or not the rearrangement is in frame (T/F) |
| stop_codon | Whether or not the rearrangement contains a stop codon (T/F) |
| conserved_c | Whether or not the rearrangement contains a detectable conserved cysteine (T/F) |
| conserved_f | Whether or not the rearrangement contains a detectable conserved phenylalanine or equivalent (T/F) |
| legacy_v_call | What older versions of Decombinator (i.e. <= v3) referred to this V gene as: e.g. 'TRBV12-3,TRBV12-4' in v4 was previously referred to just as 'TRBV12-4' | 
| legacy_j_call | What older versions of Decombinator (<= 3) referred to this J gene as | 
| v_alleles | List of V gene alleles covered by this rearrangment's tag | 
| j_alleles | List of J gene alleles covered by this rearrangment's tag | 
| v_gene_functionality | [IMGT predicted functionality](http://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html#P1-2) of V gene (or genes) used in this rearrangement (F/ORF/P, comma delimited) | 
| j_gene_functionality | [IMGT predicted functionality](http://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html#P1-2) of J gene (or genes) used in this rearrangement (F/ORF/P, comma delimited) | 
| sequence_alignment | Format required field - left blank |
| germline_alignment | Format required field - left blank |
| v_cigar | Format required field - left blank |
| d_cigar | Format required field - left blank |
| j_cigar | Format required field - left blank |
| av_UMI_cluster_size | The average UMI count for this particular sequence |

You can also use the 'nonproductivefilter' flag  (`-npf`) to suppress the output of non-productive rearrangements. 

As the hypervariable region and the primary site of antigenic contact, the CDR3 is almost certainly going to be the region of most interest for many analyses. By convention, the [CDR3 junction is defined as](http://dx.doi.org/10.1016/S0145-305X(02)00039-3) running from the position of the second conserved cysteine encoded in the 3' of the V gene to the phenylalanine in the conserved 'FGXG' motif in the J gene. However, some genes use non-canonical residues/motifs, and the position of these motifs varies.

In looking for CDR3s, we also find 'non-productive' reads, i.e. those that don't appear to be able to make productive, working TCRs. This is determined based on the presence of stop codons, being out of frame, or lacking appropriate CDR3 motifs. 

The process occurs like so:

```
# Starting with a Decombinator index
43, 5, 1, 7, AGGCAGGGATC

# Used to construct whole nucleotide sequences, using the germline FASTAs as references
GATACTGGAGTCTCCCAGAACCCCAGACACAAGATCACAAAGAGGGGACAGAATGTAACTTTCAGGTGTGATCCAATTTCTGAACACAACCGCCTTTATTGGTACCGACAGACCCTGGGGCAGGGCCCAGAGTTTCTGACTTACTTCCAGAATGAAGCTCAACTAGAAAAATCAAGGCTGCTCAGTGATCGGTTCTCTGCAGAGAGGCCTAAGGGATCTTTCTCCACCTTGGAGATCCAGCGCACAGAGCAGGGGGACTCGGCCATGTATCTCTGTGCCAGCAGCTTAGAGGCAGGGATCAATTCACCCCTCCACTTTGGGAATGGGACCAGGCTCACTGTGACAG

# This is then translated into protein sequence
DTGVSQNPRHKITKRGQNVTFRCDPISEHNRLYWYRQTLGQGPEFLTYFQNEAQLEKSRLLSDRFSAERPKGSFSTLEIQRTEQGDSAMYLCASSLEAGINSPLHFGNGTRLTVT

# The CDR3 sequence is then extracted based on the conserved C and FGXG motifs (as stored in the .translate supplementary files)
CASSLEAGINSPLHF
```

In order to do so, a third kind of supplementary data file is used, .translate files, which provide the additional information required for CDR3 extraction for each gene type (a/b/g/d, V/J). They are stored in the [TCR tag repository](https://github.com/innate2adaptive/Decombinator-Tags-FASTAs) and meet the same naming conventions as the tag and FASTA files and consist of four comma-delimited fields, detailing:
* Gene name
* Conserved motif position (whether C or FGXG)
* Conserved motif sequence (to account for the non-canonical)
* IMGT-defined gene functionality (F/ORF/P)

  
</details>

## Advanced Usage

The `decombinator` pipeline is often run over many GB of data over many hours. Under these circumstances, it can be more practical to use a cluster of workstations rather than run the pipeline locally. This section provides instructions for setting up the pipeline to run on clusters managed by University College London, but please take inspiration from our methods and adapt them to your local compute resources.

**Using the UCL CS HPC is recommended for Innate2Adaptive group members.**

See the [Decombinator-Tools](https://github.com/innate2adaptive/Decombinator-Tools/tree/master/jobs/cshpc) repository for scripts and walkthrough for group members.

To setup a similar environment for your own group, install `decombinator` into an environment, download `Decombinator-Tags-FASTAs` to an accessible location, and adapt the scripts in the above repository to be compatible with your HPC service's job engine and directory structure.

## Development Environment

```shell
git clone https://github.com/innate2adaptive/Decombinator.git
```

Install the required dependencies:

```shell
pip install .[dev]
```

Tests can be then run via:

```shell
pytest
```

To run tests offline, the `Decombinator-Tags-FASTAs` git submodule must be initialised via:

```shell
git submodule init
git submodule update
```

## How to collaborate

Please see our [code of conduct](./CODE_OF_CONDUCT.md) before contributing.

To report a bug or request support please post an issue [here](https://github.com/innate2adaptive/Decombinator/issues) outlining the problem faced and including error messages and logs where possible.

To propose a pull request, please create an issue first to discuss your proposed changes. We use [Black](https://github.com/psf/black) formatting, with a line length of 80.

For any support, commerical inquiries, or other requests, contact us [here](m.cowley@ucl.ac.uk).

## Related reading

The published history of the development of the pipeline is covered in the following publications:

* [Thomas et al (2013), Bioinformatics: *Decombinator: a tool for fast, efficient gene assignment in T-cell receptor sequences using a finite state machine*](http://dx.doi.org/10.1093/bioinformatics/btt004)
* [Oakes et al (2017), Frontiers in Immunology: *Quantitative Characterization of the T Cell Receptor Repertoire of Naïve and Memory Subsets Using an Integrated Experimental and Computational Pipeline Which Is Robust, Economical, and Versatile*](https://doi.org/10.3389/fimmu.2017.01267)
* [Uddin et al (2019), Cancer Immunosurveillance: *An Economical, Quantitative, and Robust Protocol for High-Throughput T Cell Receptor Sequencing from Tumor or Blood*](http:/dx.doi.org/10.1007/978-1-4939-8885-3_2)
* [Uddin et al (2019), Methods in Enzymology: *Quantitative analysis of the T cell receptor repertoire*](https://doi.org/10.1016/bs.mie.2019.05.054)
* [Peacock et al (2020), Bioinformatics: *Decombinator V4: an improved AIRR compliant-software package for T-cell receptor sequence annotation*](https://doi.org/10.1093/bioinformatics/btaa758)
