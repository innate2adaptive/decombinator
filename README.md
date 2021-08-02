# innate2adaptive / Decombinator 
## v4.2.0

##### Innate2Adaptive lab @ University College London, 2020
##### Written by James M. Heather, Tahel Ronel, Thomas Peacock, Niclas Thomas and Benny Chain, with help from Katharine Best, Theres Oakes and Mazlina Ismail.

---

Decombinator is a fast and efficient tool for the analysis of T-cell receptor (TCR) repertoire sequences produced by deep sequencing. The current version (v4) improves upon the previous through:
* improving the UMI based error-correction
* outputting results in a AIRR-seq Community compatible format
* increasing the dictionary of 8bp SP2 indexes from 12 to 88

## For the original Decombinator paper, please see [Thomas *et al*, Bioinformatics, 2013](http://dx.doi.org/10.1093/bioinformatics/btt004).

<h3 id="top">Navigation:</h3>

* [Installation](#installation)
* [Test data](#testdata)
* [Demultiplexing](#demultiplexing)
* [Decombining](#decombinator)
* [Collapsing](#collapsing)
* [CDR3 extraction](#cdr3)
* [Supplementary Scripts](#supplementary-scripts)
* [General usage notes](#generalusage)
* [Calling Decombinator from other scripts](#calldcr)

---

TCR repertoire sequencing (TCRseq) offers a powerful means to investigate biological samples to see what sequences are present in what distribution. However current high-throughput sequencing (HTS) technologies can produce large amounts of raw data, which presents a computational burden to analysis. Such raw data will also unavoidably contain errors relative to the original input molecules, which could confound the results of any experiment.

Decombinator addresses the problem of speed by employing a rapid and [highly efficient string matching algorithm](https://figshare.com/articles/Aho_Corasick_String_Matching_Video/771968) to search the FASTQ files produced by HTS machines for rearranged TCR sequence. The central algorithm searches for 'tag' sequences, the presence of which uniquely indicates the inclusion of particular V or J genes in a recombination. If V and J tags are found, Decombinator can then deduce where the ends of the germline V and J gene sections are (i.e. how much nucleotide removal occurred during V(D)J recombination), and what nucleotide sequence (the 'insert sequence') remains between the two. These five pieces of information - the V and J genes used, how many deletions each had and the insert sequence - contain all of the information required to reconstruct the whole TCR nucleotide sequence, in a more readily stored and analysed way. Decombinator therefore rapidly searches through FASTQ files and outputs these five fields into comma delimited output files, one five-part classifier per line.

Version 4 and higher of the Decombinator suite of scripts are written in **Python v3.7**. Older versions are written in Python v2.7. The default parameters are set to analyse data as produced by the ligation-mediated 5' RACE TCR amplification pipeline. The pipeline consists of four scripts, which are applied sequentially to the output of the previous, starting with TCR-containing FASTQ files (produced using the barcoding 5' RACE protocol):

1. Raw FASTQ files are demultiplexed to individual samples
2. Sample specific files are then searched for rearranged TCRs with Decombinator
3. Decombined data is error-corrected ('collapsed') with reference to the random barcode sequences added prior to amplification
4. Error-corrected data is translated and CDR3 sequences are extracted, in addition to the IMGT V and J gene names, the full-length variable region sequence (and the 'collapsed' sample abundance if applicable).

Very large data containing many samples, such as from Illumina NextSeq machines, often require an additional step prior to demultiplexing, to convert the raw .bcl files output by the sequencer into the three FASTQ read files (usually using `bcl2fastq`).

---

<h1 id="installation">Installation and getting started</h1>

*Please note that some familiarity with Python and command-line script usage is assumed, but everything should be easily managed by those without strong bioinformatic or computational expertise.*

### Virtual Environments

We strongly recommend running the Decombinator pipeline within a virtual environment to avoid potential package dependency clash with other installed python projects. [Conda](https://conda.io/projects/conda/en/latest/index.html) is an easy to use tool to set up separate virtual environments for your projects. This section gives some brief instructions for running Decombinator using Conda. If you are familiar with working with these tools, proceed to the [Required modules](#required-modules) section.

1. First follow the instructions provided [here](https://docs.conda.io/en/latest/miniconda.html) to install Miniconda. (Alternatively use [Anaconda](https://docs.anaconda.com/anaconda/install/) - Minconda is smaller conda distribution that is faster to install, and runs with no loss of functionality for Decombinator.)

2. Once downloaded, you can verify conda has installed correctly using the following command in your bash environment:
    ```bash
    conda info
    ```

3. If building a virtual environment for the first time, create a python v3.7 virtual environment for Decombinator using the following command:
    ```bash
    conda create --name dcrpy3 python=3.7 
    ```
    You can replace `dcrpy3` with your preferred name for the environment. **Note:** You only need to run this command once for the initial environment creation. If you already have an environment created, skip ahead to step 4.

4. To see a list of your existing conda environments, run the command:
    ```bash
    conda env list
    ```
5. Activate your environment using:
    ```bash
    source activate dcrpy3
    ```
6. Check which packages and package versions you have installed in your environment using:
    ```
    conda list
    ```
7. That's it! Now you have your environment set up, you can proceed to installing the [non-standard packages](#required-modules) required for running Decombinator.

8. To deactivate your environment, use:
    ```bash
    conda deactivate
    ```
    **Note:** this command will not destroy or delete your environment or installed packages.

### Required modules

Python 3.7 or above is required to run this pipeline, along with the following non-standard modules:

* acora (>= 2.2)
* biopython (>= 1.75)
* networkx (>= 2.5)
* polyleven (>= 0.5)
* python-Levenshtein (>= 0.12.0)
* regex (>= 2020.7.14)

**Note:** we aim to keep Decombinator up to date with the latest versions of all packages. The versions given above represent the most recently tested versions of the non-standard packages required for Decombinator.

These modules can be installed via pip (although most will likely appear in other package managers). Pip is a standard package that is automatically installed as part of Anaconda or Minconda. Install the non-standard packages by running the following command:
```bash
pip install acora>=2.2 biopython>=1.75 networkx>=2.5 polyleven>=0.5 python-levenshtein>=0.12.0 regex>=2020.7.14
```
If users are unable to install Python and the required modules, a Python installation complete with necessary packages has been bundled into a Docker container, which should be runnable on most setups. The appropriate image is located [on Dockerhub, under the name 'dcrpython'](https://hub.docker.com/r/decombinator/dcrpython/). (Please note that this is not yet updated for v4)

If you are using Windows you may need to install VS Buildtools in order to install some packages. 

### Get scripts

The Decombinator scripts are all downloadable from the [innate2adaptive/Decombinator Github repository](https://github.com/innate2adaptive/Decombinator). These can be easily downloaded in a Unix-style terminal like so:

```bash
git clone https://github.com/innate2adaptive/Decombinator.git
```

Decombinator also requires a number of additional files, which contain information regarding the V and J gene sequences, their tags, and the locations and motifs which define their CDR3 regions. By default Decombinator downloads these files from [the git repository where they are maintained](https://github.com/innate2adaptive/Decombinator-Tags-FASTAs), which obviously requires a working internet connection. In order to run Decombinator offline, these files must be downloaded to a local location, and either stored within the directory where you wish to run Decombinator or specified using the appropriate command line flag. These additional files can also be downloaded via `git clone`:

```bash
git clone https://github.com/innate2adaptive/Decombinator-Tags-FASTAs.git
```

The current version of Decombinator has tag sets for the analysis of alpha/beta and gamma/delta TCR repertoires from both mouse and man. In addition to the original tag set developed for the 2013 Thomas *et al* Bioinformatics paper, an 'extended' tag set has been developed for human alpha/beta repertoires, which covers 'non-functional' V and J genes (ORFs and pseudogenes, according to IMGT nomenclature) as well as just the functional.

### Running Decombinator on a Cluster

The Decombinator pipeline is often run over many GB of data over many hours. Under these circumstances, it can be more practical to use a cluster of workstations rather than run the pipeline locally. This section provides instructions for setting up the pipeline to run on clusters managed by University College London (Legion/Myriad).

* You will first need an account to run jobs on the cluster. Follow the instructions [here](https://www.rc.ucl.ac.uk/docs/Clusters/Myriad/).
* Once you have logged into the cluster, you should assess which python versions are currently installed:
    ```bash
    module avail python
    ```
* Look for the option including `miniconda3` in the output of the previous command. Load this module using the following command:
    ```bash
    module load python/miniconda3/4.5.11
    ```
    This will load Miniconda, and allow you to create a virtual environment for Decombinator in your local space on the cluster.
* Follow the instructions given in the [Virtual Environments](#virtual-environments) section.
Once your environment is activated, install the required non-standard packages as detailed in the [Required modules](#required-modules) section.
* Install Decombinator using git as described in the [Get Scripts](#get-scripts) section.
* To run Decombinator on the cluster, you should familiarise yourself with the process of writing, submitting and monitoring job scripts. Guidance for new users is provided [here](https://www.rc.ucl.ac.uk/docs/New_Users/).
* **Note:** you will need to include the `module load python/miniconda3/4.5.11` command and `source activate venvname` in your job scripts before calling the Decombinator scripts.
* An example job script for running the Decombinator Test Data on a cluster can be found in the `recipes/jobscripts` directory in the [Decombinator-Tools](https://github.com/innate2adaptive/Decombinator-Tools) repository. This template script can be modified to suit your own data and requirements.

### General notes

* All scripts accept gzipped or uncompressed files as input
* By default, all output files will be gz compressed
 * This can be suppressed by using the "don't zip" flag: `-dz`
* All files also output a summary log file
 * We strongly recommend you familiarise yourself with them.
* All options for a given script can be accessed by using the help flag `-h`

<sub>[↑Top](#top)</sub>

---

<h1 id="testdata">Test data</h1>

##### A small test data set of FASTQ reads produced using the Innate2Adaptive TCR amplification protocol, along with a basic set of instructions, are available [here](https://github.com/innate2adaptive/Decombinator-Test-Data).

```bash
git clone https://github.com/innate2adaptive/Decombinator-Test-Data.git
```

<sub>[↑Top](#top)</sub>

---

<h1 id="demultiplexing">Demultiplexor.py</h1>

## Demultiplexor : Demultiplexing libraries from Nextseq or Novoseq to get sample-specific, barcoded V(D)J containing reads
## Allows fuzzy demultiplexing i.e. allows a specified number of mismatches in the index sequence

## Background

*This demultiplexing step is designed specifically to make use of the random barcode sequences introduced during the Chain lab's wet lab TCR amplification protocol. While it should be fairly straightforward to adapt this to other UMI-based approaches, that will require some light modification of the scripts. Users wanting to apply Decombinator to demultiplexed, non-barcoded data can skip this step.*
* NOTE from V4.2 onwards, demultiplexor produces two outputs R1 and R2; and does not copy the beginning of R2 to the beginning of R1
* Instead, the barcode is extracted later by Decombinator; this allows for barcodes at teh beginning of R1 (new mutlplex protocol) 
* and R2 (traditional ligation protocol)


You need to provide the demultiplexing script:

* The location of the 3 or 4 read files
* An index csv file giving sample names and index seqeunces 
 * Sample names will be carried downstream, so use sensible identifiers 
 * Including the chain (e.g. 'alpha') will allow auto-detection in subsequent scripts (if *only one* chain is used per file)
 * Do not use space  or '.' characters
* See example index file indexfile_test.csv in test data for the correct format

An example command for the test data set looks like this:

```bash
python Demultiplexor.py -r1 read1_test.fq.gz -r2 read2_test.fq.gz -i1 index1_test.fq.gz -i2 index2_test.fq.gz -ix indexfile_test.csv
```

If your read files are demultiplexed by the machine, using just the SP2 indexes alone. Single read files can be produced using bash to output all appropriate reads into one file, e.g.:

```bash
# Linux machines
zcat \*R1\* | gzip > R1.fq.gz
# Macs
gunzip -c \*R1\* | gzip > R1.fq.gz
```

### INPUT 

 Requires the command line input of  3 or 4 file names, giving the two Illumina seqeunce reads, plus the one or two index reads.
 Files may be uncompressed, or gzipped (and be named accordingly, e.g. File.fastq.gz)
 A additonal optional comma delimited file detailing sample index specifics is strongly recommended, allowing 
 production of correctly named files
 File must give the following details, one sample (or index combination) per line, with no empty lines:
     Sample name, SP1/R1 index (I), SP2/R2 index (L):
     e.g.: AlphaSample,1,11
     If you include one and only one chain description (i.e. alpha, beta, gamma or delta) into your sample name, you need not set chain in Decombinator
    
** NOTE 
V4.2 simply takes the two Illumina seqeunce reads, demultiplexes and outputs two reads per sample (annotated _R1, and _R2) 
plus two undetermined reads.  In contrast to earlier versions, note that Demultiplexor does not attach a barcode from 
read 2 to read 1. This is done direcly in Decombinator.

** Other optional flags:
  
   -s/--supresssummary: Supress the production of a summary file containing details of the run into a 'Logs' directory. 
  
   -a/--outputall: Output the results of all possible index combinations currently used in protocol
*   e.g. Useful in finding potential cross-contaminating or incorrectly indexed samples
*   NB - This option can be run even if an index list is provided (although only those provided by the index list will be named)
  
   -t/--threshold: Specifies the threshold by which indexes can be clustered by fuzzy string matching, allowing for sequencing errors
*     Default = 2. Setting to zero turns off fuzzy matching, i.e. only allowing exact string matching
  
  -dz/--dontgzip: Suppress the automatic compression of output demultiplexed FASTQ files with gzip. 
*     Using this flag makes the script execute faster, but data will require more storage space. 
    
   -dc/--dontcount: Suppress the whether or not to show the running line count, every 100,000 reads. 
*     Helps in monitoring progress of large batches. 
    
   -fz/--fuzzylist: Output a list of FASTQ IDs of reads which are demultiplexed using fuzzy (i.e. non-exact) index matching, within the specified threshold.
*     Default = False, but can be useful to investigate suspected cases of poor quality index reads or clashing sequences.

  -ex/--extension: Allows users to specify the file extension of the demultiplexed FASTQ files produced.

  -cl/--compresslevel: Allows user to specify the speed of gzip compression of output files as an integer from 1 to 9. 
*     1 is the fastest but offers least compression, 9 is the slowest and offers the most compression. Default for this program is 4. 

* To see all options, run: python Demultiplexor.py -h



## OUTPUT 

    
** Versions up to 4.2. 
A fastq file will be produced for each sample listed in the index file, in the modified format, containing all reads that matched
 So we go from:        R1 - [6s|X1|----J(D)V-----]  
                       R2 - [X2]
                       R3 - [8s|N1|8s|N2|2s|-----5'UTR-----]
 To: ========>         out- [8s|N1|8s|N2|2s|X1|X2|----J(D)V-----]         
 Where X = hexamer index sequence, N = random barcode sequence, and Ys = spacer sequences of length Y
   The 8s sequences can be used downstream to identify the presence and location of random barcode N sequences
   2s is also kept to allow for the possibility finding N sequences produced from slipped reads

** Version 4.2 Produces 2 outputs per paired index, which are simply R1 and R2 Illumina reads. 
* NOTE No barcode manipulation is carried out any longer



Addition of new index sequences will currently require some slight modification of the code, although if people requested the use of a more easily edited external index file that could be incorporated in the next update.


<sub>[↑Top](#top)</sub>

---

<h1 id="decombinator">Decombinator.py</h1>

## Decombinator : identifying rearranged TCRs and outputting their 5-part classifiers, together with unique molecular identifier 

This script performs the key functions of the pipeline, as it searches through demultiplexed reads for rearranged TCR sequences. It looks for short 'tag' sequences (using Aho-Corasick string matching): the presence of a tag uniquely identifies a particular V or J gene. If it finds both a V and a J tag (and the read passes various filters), it assigns the read as recombined, and outputs a five-part Decombinator index (or 'DCR'), which uniquely represents a given TCR rearrangement.

All DCR-containing output files are comma-delimited, with the fields of that five part classifier containing, in order:
* The V index (which V gene was used)
* The J index
* Number of V deletions (relative to germ line)
* Number of J deletions
* Insert sequence (the nucleotide sequence between end of deleted V and J)

The V and J indices are arbitrary numbers based on the order of the tag sequences in the relevant tag file (using Python indexing, which starts at 0 rather than 1). Also note that the number of V and J deletions actually just represents how many bases have been removed from the end of that particular germline gene (as given in the germline FASTA files in the additional file repo); it is entirely possible that more bases where actually deleted, and just that the same bases have been re-added.

Additionally there are low frequencies of (predominantly alpha chain) recombinations where there is no detectable insertion, and where the nucleotides at the junction between the germline V and J genes could have derived from either. In such circumstances the nucleotides will arbitrarily be deemed to have derived from the V gene, and thus count towards deletions from the J, however it is impossible to know which gene originally contributed these residues.

Various additional fields may follow the five part classifier, but the DCR will always occupy the first five positions. An example identifier, from a human alpha chain file, might look like this:

```bash
1, 22, 9, 0, CTCTA
```

Which corresponds to a rearrangement between TRAV1-2 (V index **1**, with **9** nucleotides deleted) and TRAJ33 (J index **22**, with **0** deletions), with an insert sequence (i.e. non-templated additions to the V and/or the J gene) of '**CTCTA**'. For beta chains, the insert sequence will contain any residual TRBD nucleotides, although as these genes are very short, homologous and typically highly 'nibbled', they are often impossible to differentiate.

```bash
python Decombinator.py -fq AlphaSample1.fq.gz -sp mouse -c a
```
** Version 4.2 introduces some changes to work with Demultiplexor 4.01 and to work for the new multiplex protocol. It can look for a 
** barcode in either read 1 (multiplex protocol) or read 2 (ligation protocol). This is controled by a new required flag -bc_read 
** which must be R1 or R2. The bc_length can also be set - the default is 42. An example command for the multiplex protocol may look like this:

```bash
python Decombinator.py -fq BetaSample1.fq.gz -c b -br R1 -bl 22
```


## INPUT 

 As with entire pipeline, Decombintator is run using command line arguments to provide user parameters
   All arguments can be read by viewing the help data, by running python Decombintator.py -h

** The two required parameters are 
 1. -fq/--fastq which identify FASTQ reads produced by Demultiplexor.py (unzipped or gzipped).
 2 -br/--bc_read which determines whether teh barcode is obtained from the beginning of R2 as in the standard 5'RACE ligation protocol, using the M13-I8-6N-I8-6N or the older SP2-I8-6N-I8-6N oligonucleotide; or from the beginning of R1 as in the new Vbeta multiplex protocol. The length of sequence containing the barcode can also be determined using -bl/--bc_length. The default is 42; but for the multiplex, 22 is enough. 

 The TCR chain locus to look for can be explicitly specified using the -c flag 
 Users can use their choice of chain identifiers from this list (case insensitive): a/b/g/d/alpha/beta/gamma/delta/TRA/TRB/TRG/TRD/TCRA/TCRB/TCRG/TCRD
 If no chain is provided (or if users which to minimise input arguments), script can infer chain from the FASTQ filename
 I.e. "alpha_sample.fq" would be searched for alpha chain recombinations
** NB: This autodetection only works if there is only ONE TCR locus present in name (which must be spelt out in full)

** Other optional flags:
  
   -s/--supresssummary: Supress the production of a summary file containing details of the run into a 'Logs' directory. 
      
   -dz/--dontgzip: Suppress the automatic compression of output demultiplexed FASTQ files with gzip. 
*  Using this flag makes the script execute faster, but data will require more storage space. 
    
   -dc/--dontcount: Suppress the whether or not to show the running line count, every 100,000 reads. 
*  Helps in monitoring progress of large batches.
  
   -dk/--dontcheck: Suppress the FASTQ sanity check. 
 * Strongly recommended to leave alone: sanity check inspects first FASTQ read for basic FASTQ parameters.
  
   -pf/--prefix: Allows users to specify the prefix of the Decombinator TCR index files produced. Default = 'dcr_'
  
   -ex/--extension: Allows users to specify the file extension of the Decombinator TCR index files produced. Default = '.n12'

  -or/--orientation: Allows users to specify which DNA orientations to check for TCR reads. Default = reverse only, as that's what the protocol produces.
 * This will likely need to be changed for analysing data produced by protocols other than our own.

   -tg/--tags: Allows users to specify which tag set they wish to use. For human alpha/beta TCRs, a new 'extended' tag set is recommended, as it covers more genes.
 *   An extended tag set is only currently available for human a/b genes.

  -sp/--species: Current options are only human or mouse. Help could potentially be provided for generation of tags for different species upon request.
  
  -N/--allowNs: Provides users the option to allow 'N's (ambiguous base calls), overriding the filter that typically removes rearrangements that contain them.
 *  Users are recommended to not allow Ns, as such bases are both themselves low quality data and predict reads that are generally less trustworthy.
    
 * -ln/--lenthreshold: The length threshold which (the inter-tag region of) successful rearrangements must be under to be accepted. Default = 130.
  
   -tfdir/--tagfastadir: The path to a local copy of a folder containing the FASTA and Decombinator tag files required for offline analysis.
 *    Ordinarily such files can be downloaded on the fly, reducing local clutter. By default the script looks for the required files in the present working directory, then in a subdirectory called "Decombinator-Tags-FASTAs", then online.
 *    Files are hosted on GitHub, here: https://github.com/innate2adaptive/Decombinator-Tags-FASTAs

   -nbc/--nobarcoding: Run Decombinator without any barcoding, i.e. use the whole read. 
 * Recommended when running on data not produced using the Innate2Adaptive lab's ligation-mediated amplification protocol
  
  -bl/--bc_length : sets the length of seqeunce to be stored by Decombinator from R1 or R2 (as set by -bc_read) for further use by Collapsinator.
   

  ## OUTPUT 
  
  
  Produces a '.n12' file by default, which is a standard comma-delimited Decombinator output file with several additional fields:
  V index, J index, # V deletions, # J deletions, insert, ID, TCR sequence, TCR quality, barcode sequence, barcode quality
**  NB The TCR sequence given here is the 'inter-tag' region, i.e. the sequence between the start of the found V tag the end of the found J tag 

<sub>[↑Top](#top)</sub>

---

<h1 id="collapsing">Collapsinator.py</h1>

## Collapsing: using the random barcodes to error- and frequency-correct the repertoire

Takes the output files of Decombinator (run using the barcoding option) and performs collapsing and error correction. This version is a modified version of KB's script collapsinator_20141126.py (That was itself an improved version of the CollapseTCRs.py script used in the Heather et al HIV TCR paper (DOI: 10.3389/fimmu.2015.00644))*
**  Version 4.0.2 includes improved clustering routines measuring the similarity in both barcode and TCR sequence of TCR repertoire data
  
**  NOTE - from version 4.2 this optionally looks for barcode 6NI86N at the beginning of the read; instead of M13_6N_I8_6N_I8
  (i.e. only one spacer).
  This makes it compatible with the multiplex protocol in which the barcode is incorproated in the RT step and is found at the beginning of R1. 
 ** from version V4.2  there is a required additional command line parameter -ol (see below for allowed inputs)
The barcode sequence is contained in one of the additional fields output by `Decombinator.py` in the .n12 files, that which contains the first 42 bases of R2. As Illumina sequencing is particularly error-prone in the reverse read, and that reads can be phased (i.e. they do not always begin with the next nucleotide that follows the sequencing primer) our protocol uses known spacer sequences to border the random barcode bases, so that we can identify the actual random bases. The hexameric barcode locations (N6) are determined in reference to the two spacer sequences like so:

```
I8 (spacer) – N6 – I8 – N6 – 2 base overflow (n)
GTCGTGATNNNNNNGTCGTGATNNNNNNnn
```

The collapsing script uses the spacer sequences to identify the exact position of the barcode sequences.

The `Collapsinator.py` script performs the following procedures:
* Scrolls through each line of the input .n12 file containing DCR, barcode and sequence data.
* Removes TCR reads with forbidden errors, e.g. ambigious base calls (with user input parameters provided to modify strictness).
* Groups input reads by barcode. Reads with identical barcode and equivalent inter-tag sequence are grouped together. Equivalence is defined as the Levenshtein distance between two sequences being lower than a given threshold, weighted by the lengths of the compared sequences. Reads with identical barcodes but non-equivalent sequences are grouped separately.
* Each group is assigned the most common inter-tag sequence/DCR combination as the 'true' TCR, as errors are likely to occur during later PCR cycles, and thus will most often be minority variants (see [Bolotin *et al.*, 2012](http://dx.doi.org/10.1002/eji.201242517)).

After this initial grouping, the script estimates the true cDNA frequency. UMIs that are both similar and are associated to a similar TCR are likely to be amplified from the same original DNA molecule, and to differ only due to PCR or sequencing error. Consequently, groups with similar barcodes and sequences are then clustered via the following procedure:
* The barcode of each group is compared to the barcode of every other group.
* The expected distribution of distances between UMIs can be modelled as a binomial distribution. Experimentation with simulated datasets found the best threshold for allowing two barcodes to be considered equivalent is when they have a Levenshtein distance of less than 3; a value of 2 is set by default. This can be modified through the user input parameter ```-bc```.
* Groups with barcodes that meet this threshold criteria have their inter-tag sequences compared. Those with equivalent sequences are clustered together. Sequence equivalence is here taken to mean that the two sequences have a Levenshtein distance less than or equal to 10% of the length of the shorter of the two sequences. This percentage can be modified through the user input parameter ```-lv```.
* Upon this merging of groups, the most common inter-tag sequence of the cluster is reassessed and taken as the 'true' TCR. 

Finally, the clusters are collapsed to give the abundance of each TCR in the biological sample.
* A TCR abundance count is calculated for each TCR by counting the number of clusters that have the same sequence but different barcodes (thus representing the same rearrangement originating from multiple input DNA molecules).
* An average UMI count is calculated for each TCR by summing the number of members in each cluster associated with the TCR sequence, and dividing by the number of those clusters. This gives a measure that can be used to estimate the robustness of the data for that particular sequence.

Collapsinator outputs 7 fields: the 5-part DCR identifier, the corrected abundance of that TCR in the sample, and the average UMI count for that TCR

Collapsing occurs in a chain-blind manner, and so only the decombined '.n12' file is required, without any chain designation, with the only required parameter being the infile:

```bash
python Collapsinator.py -in dcr_AlphaSample1.n12
```

From Version 4.2 onwards, the oligo flag needs to be included in the command. For example, to run Collapsinator on samples from the multiplex protocol

```bash
python Collapsinator.py -in dcr_BetaSample1.n12 -ol I8_single
```

A number of the filters and thresholds can be altered using different command line flags. In particular, changing the R2 barcode quality score and TCR sequence edit distance thresholds (via the `-mq` `-bm` `-aq` and `-lv` flags) are the most influential parameters. However the need for such fine tuning will likely be very protocol-specific, and is only suggested for advanced users, and with careful data validation. A histogram of the average UMI counts can be generated using the `-uh` flag.

The default file extension is '.freq'.

## Collapsinator Test Suite

A large number of unit tests for the Collapsinator error correction stage can be found in the `unittests` folder. These tests cover a range of cases to ensure the correct analysis of spacers with minor errors. Separate sets of tests are included for both the former and current experimental protocols used by the Chain Lab. These can act as a template for anyone wishing to design tests for differing protocols.

Individual tests can be run as in the following example:
```
python unittests/I8tests/I8correctlength.py
```
Additionally, a script is included to run all current tests at once:
```
python unittests/alltests.py
```
## INPUT 

  
**  Required inputs 
  -in/--infile : Defines input file.   Takes as input .n12 files produced by Decombinator (v3 or higher), assuming it has been run on suitably barcoded and demultiplexed data.
 
  -ol/--oligo  : Specifies the spacer (protocol dependent) as M13, I8, I8_single. The I8 protocol is deprecated.
  
**  Other optional flags:
  
    -s/--supresssummary: Supress the production of a summary file containing details of the run into a 'Logs' directory. 
  
    -dz/--dontgzip: Suppress the automatic compression of output demultiplexed FASTQ files with gzip. 
  
    -dc/--dontcount: Suppress the whether or not to show the running line count, every 100,000 reads. Helps in monitoring progress of large batches.

  The other optional flags are somewhat complex, and caution is advised in their alteration.

  To see all options, run: python Collapsinator.py -h

  Input files need to be in the appropriate format, consisting of:
    V index, J index, V deletions, J deletions, insert, ID, inter-tag TCR sequence, inter-tag quality, barcode sequence, barcode quality


## OUTPUT 
   
  A Decombinator index file, giving each error-corrected DCR index, and the frequency with which it appears
  in the final processed data, and an average UMI count, which can be used to estimate the robustness of the
  data for that particular sequence
<sub>[↑Top](#top)</sub>

---

<h1 id="cdr3">CDR3translator.py</h1>

## TCR translation and CDR3 extraction: turning DCR indexes into complementarity determining region 3 sequences

As the hypervariable region and the primary site of antigenic contact, the CDR3 is almost certainly going to be the region of most interest for most analyses. By convention, the [CDR3 junction is defined](http://dx.doi.org/10.1016/S0145-305X(02)00039-3) as running from the position of the second conserved cysteine encoded in the 3' of the V gene to the phenylalanine in the conserved 'FGXG' motif in the J gene. However, some genes use non-canonical residues/motifs, and the position of these motifs varies.

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

In order to do so, a third kind of supplementary data file is used, .translate files, which provide the additional information required for CDR3 extraction for each gene type (a/b/g/d, V/J). They meet the same naming conventions as the tag and FASTA files, and consist of four comma-delimited fields, detailing:
* Gene name
* Conserved motif position (whether C or FGXG)
* Conserved motif sequence (to account for the the non-canonical)
* IMGT-defined gene functionality (F/ORF/P)

`CDR3translator.py` only requires an input file, and to be told which TCR chain is being examined (although this can also be inferred from the file name, as with `Decombinator.py`. Input can be collapsed or not (e.g. .freq or .n12), and the program will retain the frequency information if present.

```bash
python CDR3translator.py -in dcr_AlphaSample1.freq
# or
python CDR3translator.py -in dcr_AnotherSample1.freq -c b
```

** NOTE As of version 4, this script now outputs a tab separated file compatible with the AIRR-seq community format, to encourage data re-use and cross-tool compatibility and comparisons. For details please see [Vander Haiden *et al.* (2018)](http://dx.doi.org/10.3389/fimmu.2018.02206) and the [AIRR community standards](https://docs.airr-community.org/). Note that this format expects certain columns to be present even if the fields are not applicable, so CDR3translator leaves these fields empty. Further fields have been added.

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

You can also use the 'nonproductivefilter' flag  (`-npf`) to suppress the output of non-productive rearrangements. 

<sub>[↑Top](#top)</sub>

---

<h1 id="supplementary-scripts">Supplementary Scripts</h1>

A number of useful scripts to supplement the Decombinator pipeline can be found in the [Decombinator-Tools](https://github.com/innate2adaptive/Decombinator-Tools) repository.

This repository includes scripts to automatically generate test data, randomly subsample existing data, automatically run the Decombinator test data, generate a histograms of average UMI cluster sizes, and to build a summary of log files.

<sub>[↑Top](#top)</sub>

---

<h1 id="generalusage">General usage notes and tips</h1>

## Being efficient

* Use the help flag (`-h`) to see all options for a given script
* Where possible, run scripts in same folder as data
    * Should work remotely, but recommended not to
* Use bash loops to iterate over many samples
    * For example, to run the CDR3 translation script on all alpha .freq files in the current directory:
`for x in *alpha*.freq.gz; do echo $x; python ligTCRtranslateCDR3.py -in $x -c a; done`
* There is also a [makefile available](https://github.com/innate2adaptive/automate-decombinator) to automate the pipeline from start to finish
* Ideally don't pool alpha and beta of one individual into same index for sequencing
    * You may wish to know the fraction of reads Decombining, which will be confounded if there are multiple samples or chains per index combination

## Be aware of the defaults

* These scripts are all set up to analyse data produced from our ligation TCR protocol, and thus the defaults reflect this
    * By far the most common libraries analysed to date are human αβ repertoires
* If running on mouse and/or gamma/delta repertoires, it may be worth altering the defaults towards the top of the code to save on having to repeatedly have to set the defaults
* If running on data not produced using our protocol, likely only the Decombinator and translation scripts will be useful 
    * The `-nbc` (non-barcoding) flag will likely need to be set (unless you perform analogous UMI-positioning to `Demultiplexor.py`)
    * If your libraries contain TCRs in the forward orientation, or in both, this will also need to be set with the `-or` flag

## Keep track of what's going on

* Using the 'dontgzip' flag ('-dz') will stop the scripts automatically compressing the output data, which will minimise processing time (but take more space)
    * It is recommend to keep gzipped options as default
    * Compressing/decompressing data slightly slows your analysis, but dramatically reduces the data storage footprint
* Don't be afraid to look at the code inside the scripts to figure out what's going on
* Check your log files, as they contain important information
* **Back up your data** (including your log files)!

## Use good naming conventions

* Keep separate alpha and beta files, named appropriately
    * Keep chain type in name
    * Decombinator will always take specified chain as the actual, but failing that will look in filename
* Don't use full stops in file names
    * Some of the scripts use them to determine file extension locations
    * Additionally discourages inadvertent data changing
* While there are default file extensions, these settings can be set via the `-ex` flag 
    * Similarly, when Decombining you can set the prefix with `-pf`
 
## Using Docker

(Please note, this is not yet updated to v4)
* Docker may be convenient for systems where installing many different packages is problematic
    * It should allow users to run Decombinator scripts without installing Python + packages (although you still have to install Docker)
    * The scripts however will run faster if you natively install Python on your OS
* If opted for, users must first [install Docker](https://docs.docker.com/engine/installation/)
* Then scripts can be run like so:

```bash
docker run -it --rm --name dcr -v "$PWD":/usr/src/myapp -w /usr/src/myapp decombinator/dcrpython python SCRIPT.py`
```

* Depending on your system, you may need to run this as a superuser, i.e. begin the command `sudo docker...`
* N.B.: the first time this is run Docker will download [the Decombinator image from Dockerhub](https://hub.docker.com/r/decombinator/dcrpython/), which will make it take longer than usual

<sup>[↑Top](#top)</sup>

---

<h1 id="calldcr">Calling Decombinator from other Python scripts</h1>

As an open source Python script, `Decombinator` (and associated functions) can be called by other scripts, allowing advanced users to incorporate our rapid TCR assignation into their own scripts, permitting the development of highly bespoke functionalities.

If you would like to call Decombinator from other scripts you need to make sure you fulfil the correct requirements upstream of actually looking for rearrangements. You will thus need to:

* Set up some 'dummy' input arguments, as there are some that Decombinator will expect
* Run the 'import_tcr_info' function, which sets up your environment with all of the required variables (by reading in the tags and FASTQs, building the tries etc)

An example script which calls Decombinator might therefore look like this:

```python
import Decombinator as dcr

# Set up dummy command line arguments 
inputargs = {"chain":"b", "species":"human", "tags":"extended", "tagfastadir":"no", "lenthreshold":130, "allowNs":False, "fastq":"blank"}

# Initalise variables
dcr.import_tcr_info(inputargs)

testseq = "CACTCTGAAGATCCAGCGCACACAGCAGGAGGACTCCGCCGTGTATCTCTGTGCCAGCAGCTTATTAGTGTTAGCGAGCTCCTACAATGAGCAGTTCTTCGGGCCAGGGACACGGCTCACCGTGCTAGAGGACCTGAAA"

# Try and Decombine this test sequence
print dcr.dcr(testseq, inputargs)
```

This script produces the following list output:
`[42, 6, 2, 0, 'TTAGTGTTAGCGAG', 22, 118]`

The first five fields of this list correspond to the standard DCR index as described above, while the last two indicate the start position of the V tag and the end position of the J tag (thus definining the 'inter-tag region').

Users may wish to therefore loop through their data that requires Decombining, and call this function on each iteration. However bear in mind that run in this manner Decombinator only checks whatever DNA orientation it's given: if you need to check both orientations then I suggest running one orientation through first, checking whether there's output, and if not then running the reverse orientation through again (which can be easily obtained using the Seq functionality of the Biopython package).

<sup>[↑Top](#top)</sup>

---

### Legacy Decombinator versions

If users wish to view previous versions of Decombinator, v2.2 is available from the old GitHub repo [uclinfectionimmunity / Decombinator](https://github.com/uclinfectionimmunity/Decombinator/). However it is not recommended that this script be used for analysis as it lacks some key error-reduction features that have been integrated into subsequent versions, and is no longer supported.

---

### Related reading

The published history of the development of the pipeline is covered in the following publications:

* [Thomas et al (2013), Bioinformatics: *Decombinator: a tool for fast, efficient gene assignment in T-cell receptor sequences using a finite state machine*](http://dx.doi.org/10.1093/bioinformatics/btt004)
* [Oakes et al (2017), Frontiers in Immunology: *Quantitative Characterization of the T Cell Receptor Repertoire of Naïve and Memory Subsets Using an Integrated Experimental and Computational Pipeline Which Is Robust, Economical, and Versatile*](https://doi.org/10.3389/fimmu.2017.01267)
* [Uddin et al (2019), Cancer Immunosurveillance: *An Economical, Quantitative, and Robust Protocol for High-Throughput T Cell Receptor Sequencing from Tumor or Blood*](http:/dx.doi.org/10.1007/978-1-4939-8885-3_2)
* [Uddin et al (2019), Methods in Enzymology: *Quantitative analysis of the T cell receptor repertoire*](https://doi.org/10.1016/bs.mie.2019.05.054)
* [Peacock et al (2020), Bioinformatics: *Decombinator V4: an improved AIRR compliant-software package for T-cell receptor sequence annotation*](https://doi.org/10.1093/bioinformatics/btaa758)

<sup>[↑Top](#top)</sup>
