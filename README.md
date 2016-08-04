# innate2adaptive / Decombinator 
## v3.1
#### This version written by James M. Heather, Innate2Adaptive lab under Benny Chain @ University College London, 2016

Decombinator is a fast and efficient tool for the analysis of T-cell receptor (TCR) repertoire sequences produced by deep sequencing. The current version (v3) improves upon the previous through:
* incorporation of error-correction strategies using unique molecular index (UMI) barcoding techniques
* optimisation for additional speed and greater accuracy
* usage of extended tag sets of V and J genes allowing detection of genes of all functionalities
* greater command line input, log file production and code commenting, providing greater flexibility, repeatability and ease of modification

## For the original Decombinator paper, please see [Thomas *et al*, Bioinformatics, 2013](http://dx.doi.org/10.1093/bioinformatics/btt004).

TCR repertoire sequencing (TCRseq) offers a powerful means to investigate biological samples to see what sequences are present in what distribution. However current high-throughput sequencing (HTS) technologies can produce large amounts of raw data, which presents a computational burden to analysis. This such raw data will also unavoidably contain errors relative to the original input molecules, which could confound the results of any experiment.

Decombinator addresses the problem of speed by employing a rapid and [highly efficient string matching algorithm](https://figshare.com/articles/Aho_Corasick_String_Matching_Video/771968) to search the FASTQ files produced by HTS machines for rearranged TCR sequence. The central algorithm searches for 'tag' sequences, the presence of which uniquely indicates the inclusion of particular V or J genes in a recombination. If V and J tags are found, Decombinator can then deduce where the ends of the germline V and J gene sections are (i.e. how much nucleotide removal occurred during V(D)J recombination), and what nucleotide sequence (the 'insert sequence') remains between the two. These five pieces of information - the V and J genes used, how many deletions each had and the insert sequene - contain all of the information required to reconstruct the whole TCR nucleotide sequence, in a more readily stored and analysed way. Decombinator therefore rapidly searches through FASTQ files and outputs these five fields into comma delimited output files, one five-part classifier per line.

The Decombinator suite of scripts are written in **Python v2.7**, and the default parameters are set to analyse data as produced by the ligation-mediated 5' RACE TCR amplification pipeline. The pipeline consists of four scripts, which are applied sequentially to the output of the previous starting with TCR-containing FASTQ files (produced using the barcoding 5' RACE protocol):

1. Raw FASTQ files demultiplexed to individual samples
2. Sample specific files searched for rearranged TCRs with Decombinator
3. Decombined data is error-corrected ('collapsed') with reference to the random barcode sequences added prior to amplification
4. Error-corrected data is translated and CDR3 sequences extracted

# Installation and getting started

*Please note that some familiarity with Python and command-line script usage is assumed, but everything should be easily managed by those without strong bioinformatic or computational expertise.*

### Required modules

Python 2.7 is required to run this pipeline, along with the following non-standard modules:

* acora
* biopython
* regex
* python-Levenshtein

These can be installed via pip (although most will likely appear in other package managers), e.g.:
```bash
pip install python-biopython python-levenshtein regex acora
```
If users are unable to install Python and the required modules, a Python installation complete with necessary packages has been bundled into a Docker container, which should be runnable on most setups. The appropriate image is located [on Dockerhub, under the name 'dcrpython'](https://hub.docker.com/r/decombinator/dcrpython/).

### Get scripts

The Decombinator scripts are all downloadable from the [innate2adaptive/Decombinator Github repository](https://github.com/innate2adaptive/Decombinator). These can be easily downloaded in a Unix-style terminal like so:

```bash
git clone https://github.com/innate2adaptive/Decombinator.git
```

Decombinator also requires a number of additional files, which contain information regarding the V and J gene sequences, their tags, and the locations and motifs which define their CDR3 regions. By default Decombinator downloads these files from [the git repository where they are maintained](https://github.com/innate2adaptive/Decombinator-Tags-FASTAs), which obviously requires a working internet connection. In order to run Decombinator offline, these files must be downloaded to a local location, and either stored within the directory where you wish to run Decombinator or specified using the appropriate command line flag. These additional files can also be downloaded via `git clone`:

```bash
git clone https://github.com/innate2adaptive/Decombinator-Tags-FASTAs.git
```

The current version of Decombinator has tag sets for the analysis of alpha/beta and gamma/beta TCR repertoires from both mouse and man. In addition to the original tag set developed for the 2013 Thomas *et al* Bioinformatics paper, an 'extended' tag set has been developed for human alpha/beta repertoires, which covers 'non-functional' V and J genes (ORFs and pseudogenes, according to IMGT nomenclature) as well as just the functional.


# Demultiplexor.py
## Demultiplexing: getting sample-specific, barcoded V(D)J containing reads

*This demultiplexing step is designed specifically to make use of the random barcode sequences introduced during the Chain lab's wet lab TCR amplification protocol. While it should be fairly straightforward to adapt this to other UMI-based approaches, that will require some light modification of the scripts. Users wanting to apply Decombinator to demultiplexed, non-barcoded data can skip this step.*

Illumina machines produce 3 FASTQ read files:

* Read 1 
 * Contains the V(D)J sequence and a demultiplexing index (the SP1 index)
* Read 2 
 * Contains the barcode sequences, and reads into the 5' UTR
* Index 1 
 * Contains the second demultiplexing index (the typical Illumina index read, giving the SP2 index)
 
The demultiplexing script extracts the relevant sequences from each read file and combines them into reads in one file, and demultiplexes them on the two indexes to produce separate sample files containing both the V(D)J and barcode information.

Note that you may need to run `bcl2fastq` or alter a configuration file on your sequencing machine in order to obtain the I1 index read file (which is not advised unless you know what you're doing, please consult whoever runs your machines). For example, this can be achieved on a MiSeq by adding the following line to the `Reporter.exe.config`, under `<appsettings>`:

```bash
<add key=”CreateFastqForIndexReads” value=”1”/>
```

You need to provide the demultiplexing script:

* The location of the 3 read files
* An index csv file giving sample names and indexes
 * Sample names will be carried downstream, so use sensible identifiers (include chain, and no '.' characters)
* Give numeric indexes, SP1 end first, then SP2 end
 * e.g. `'AlphaSample1,5,9'` would put all reads amplified using SP1-6N-I-5-aRC1 and P7-L-9 in one file

An example command might look like this:
```bash
python Demultiplexor.py -r1 R1.fq.gz -r2 R2.fq.gz -i1 I1.fq.gz -ix IndexFile.ndx
```

It's also likely that your read files are already demultiplexed by the machine, using just the SP2 indexes alone. Single read files can be produced using bash to output all appropriate reads into one file, e.g.:

```bash
# Linux machines
zcat \*R1\* | gzip > R1.fq.gz
# Macs
gunzip -c \*R1\* | gzip > R1.fq.gz
```

If you suspect your indexing might be wrong you can use the `'outputall'` flag (`-a True`). This outputs demultiplexed files using all possible combination of indexes used in our protocol, and can be used to identify the actual index combinations contained in a run and locate possible cross-contamination. Note that this is only recommended for troubleshooting and not for standard use as it will decrease the number of successfully demultiplexed reads per samples due to fuzzy index matching. 

Addition of new index sequences will currently require some slight modification of the code, although if people requested the use of a more easily edited external index file that could be incorporated in the next update.

# Decombinator.py
## Decombining: identifying rearranged TCRs and outputting their 5-part classifiers

Searches through demultiplexed reads for rearranged TCR sequences
Looks for short 'tag' sequences, the presence of which uniquely identifies a particular V or J gene
Uses Aho-Corasick string matching
If it finds both a V and a J tag (and the read passes various filters), it assigns the read as recombined, and outputs a five-part Decombinator index (or 'DCR')

5-part classifier that uniquely represents a given TCR rearrangement
Describes:
V index (which V gene was used)
J index
V deletions (number of deletions from germ line)
J deletions
Insert sequence (nucleotide sequence between end of deleted V and J)
Additional fields output for collapsing (see later)

Need to provide Decombinator a fastq file
Produced with our protocol and processed by ligTCRdemultiplex.py
Can also specify species (defaults to human) and a TCR chain (a/b/g/d)
If files were named appropriately during demultiplexing, Decombinator can detect chain without needing to specify it
python ligTCRDCR.py -fq AlphaSample1.fq.gz -sp mouse -c a
By default assigns a 'dcr_' prefix and '.n12' file extension
Aid in command line data handling, but can be altered

Decombinator uses tags to identify different V/J genes, and can therefore only detect those genes that went into the tag set
Both human and mouse have an 'original' tag set, which contains all of the prototypical alleles for each 'functional' TCR gene
As defined by IMGT
There is also an 'extended' tag set for humans, which includes tags for the prototypical alleles of all genes, regardless of predicted functionality
Recommended, due to increased specificity

Supplementary files required by Decombinator can be downloaded by the script from GitHub as it runs
If running offline, you need to download the files to the directory where you run Decombinator
Files are named by the following convention:
[species]_[tag set]_[gene].[file type]

Note the TCR gene naming convention:
TRAV1-2*01
TR = TCR / A = alpha / V = variable / 1 = family / -2 = subfamily / *01 = allele

# Collapsing


# CDR3 extractions



REMEMBER TO LINK TO MAZLINA'S MAKEFILE

NEED TO INCLUDE BIT ON USING DOCKER!!
. This should allow users to run Decombinator scripts without installing Python + packages
If opted for, users must first [install Docker](https://docs.docker.com/engine/installation/)
Then scripts can be run like so:
docker run -it --rm --name dcr -v "$PWD":/usr/src/myapp -w /usr/src/myapp decombinator/dcrpython python SCRIPT.py
You may need to write 'sudo' before 'docker'
N.B.: the first time this is run Docker will download [the Decombinator image from Dockerhub](https://hub.docker.com/r/decombinator/dcrpython/), which should have been done prior to this session





### Legacy Decombinator versions

If users wish to view previous versions of Decombinator, v2.2 is available from the old GitHub repo [uclinfectionimmunity / Decombinator](https://github.com/uclinfectionimmunity/Decombinator/). However it is not recommended that this script be used for analysis as it lacks some key error-reduction features that have been integrated to subsequent versions, and is no longer supported.




