<div align="center">

DECOMBINATOR
============

<!-- [![Latest release](https://img.shields.io/pypi/v/sceptr)](https://pypi.org/p/sceptr) -->
<!-- ![Tests](https://github.com/innate2adaptive/Decombinator/refactor/.github/workflows/package-build-test.yml/badge.svg) -->
[![License](https://img.shields.io/badge/license-MIT-blue)](https://raw.githubusercontent.com/innate2adaptive/Decombinator/master/LICENCE)
[![Code style](https://img.shields.io/badge/formatted%20with-black-black)](https://github.com/psf/black)
<!-- [![Python versions](https://img.shields.io/pypi/pyversions/pytest-docker)](https://pypi.org/project/pytest-docker/) -->
<!-- Docs -->
<!-- ### Check out the [documentation page](https://sceptr.readthedocs.io). -->

</div>

**Innate2Adaptive lab @ University College London, 2024**

**Written by (in alphabetical order) Benny Chain, James M. Heather, Katharine Best, Matthew V. Cowley, Mazlina Ismail, Niclas Thomas, Tahel Ronel, Theres Oakes, and Thomas Peacock.**

Decombinator is a fast and efficient tool for the analysis of T-cell receptor (TCR) repertoire sequences produced by deep sequencing.
It accepts TCR repertoire sequencing (TCRseq) FASTQ data and returns [AIRR](https://docs.airr-community.org/en/stable/datarep/rearrangements.html) compliant files detailing TCR V(D)J recombination and counts.

TCRseq offers a powerful means to investigate biological samples to observe sequence distributions.
However current high-throughput sequencing (HTS) technologies can produce large amounts of raw data which will also unavoidably contain errors relative to the original input molecules.
Decombinator addresses the problem of large datasets through speed - employing a rapid and [highly efficient string matching algorithm](https://figshare.com/articles/Aho_Corasick_String_Matching_Video/771968) to search the FASTQ files produced by HTS machines for rearranged TCR sequences.
The central algorithm searches for 'tag' sequences, the presence of which uniquely indicates the inclusion of particular V or J genes in a recombination.
If V and J tags are found, Decombinator can then deduce where the ends of the germline V and J gene sections are (i.e. how much nucleotide removal occurred during V(D)J recombination), and what nucleotide sequence (the 'insert sequence') remains between the two.
These five pieces of information - the V and J genes used, how many deletions each had and the insert sequence - contain all of the information required to reconstruct the whole TCR nucleotide sequence, in a more readily stored and analysed way.
The pipeline also handles FASTQ data with many reads of the same molecule due to over-amplification.
It does this by using unique molecular identifiers to "collapse" down duplicate reads.

Running decombinator is easy. All that is required is HTS read data and a few arguments specifying the chemistry of your prepared library:

```shell
decombinator pipeline -in XXXX.fq -c b -br R2 -bl 42 -ol M13
```

Installation
------------

To install decombinator and all required packages:

```shell
pip install decombinator
```

Decombinator also requires a number of additional files, which contain information regarding the V and J gene sequences, their tags, and the locations and motifs which define their CDR3 regions. By default Decombinator downloads these files from [the git repository where they are maintained](https://github.com/innate2adaptive/Decombinator-Tags-FASTAs), which obviously requires a working internet connection. In order to run Decombinator offline, these files must be downloaded to a local location, and either stored within the directory where you wish to run Decombinator or specified using the appropriate command line flag (`-tfdir`). These additional files can also be downloaded via `git clone`:

```bash
git clone https://github.com/innate2adaptive/Decombinator-Tags-FASTAs.git
```

The current version of Decombinator has tag sets for the analysis of alpha/beta and gamma/delta TCR repertoires from both mouse and man. In addition to the original tag set developed for the 2013 Thomas *et al* Bioinformatics paper, an 'extended' tag set has been developed for human alpha/beta repertoires, which covers 'non-functional' V and J genes (ORFs and pseudogenes, according to IMGT nomenclature) as well as just the functional.

Basic Usage
-----------

Decombinator can be run in 4 different modes depending on your use case. Firstly, `pipeline` mode can be used via:

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
| `-br BC_READ`, `--bc_read BC_READ`   | Which read has bar code (R1,R2)                                                                                                                                |
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
| `-ol OLIGO`, `--oligo OLIGO`         | Choose experimental oligo for correct identification of spacers ["M13", "I8","I8_single] (default: M13)                                                         |
| `-wc`, `--writeclusters`             | Write cluster data to separate cluster files                                                                                                                   |
| `-uh`, `--UMIhistogram`              | Creates histogram of average UMI cluster sizes                                                                                                                 |
| `-npf`, `--nonproductivefilter`      | Filter out non-productive reads from the output                                                                                                                |


</details>

In this mode, all three main components of the pipeline are applied to the data sequentially: `decombine`, `collapse`, and `translate`.
For full detail of each of these components please read the [documentation](tbd).

If instead, you wish to run just one component of the pipeline, these are all directly accessible via their respective sub-parser:

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
| `-br BC_READ`, `--bc_read BC_READ`   | Which read has bar code (R1,R2)                                                                                                                                |
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
| `-ol OLIGO`, `--oligo OLIGO`         | Choose experimental oligo for correct identification of spacers ["M13", "I8","I8_single] (default: M13)                                                         |
| `-wc`, `--writeclusters`             | Write cluster data to separate cluster files                                                                                                                   |
| `-uh`, `--UMIhistogram`              | Creates histogram of average UMI cluster sizes                                                                                                                 |


</details>

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

Advanced Usage
--------------

The `decombinator` pipeline is often run over many GB of data over many hours. Under these circumstances, it can be more practical to use a cluster of workstations rather than run the pipeline locally. This section provides instructions for setting up the pipeline to run on clusters managed by University College London, but please take inspiration from our methods and adapt them to your local compute resources.

**Using the UCL CS HPC is recommended for Innate2Adaptive group members.**

See the [Decombinator-Tools](https://github.com/innate2adaptive/Decombinator-Tools/tree/master/jobs/v4.3_cs_cluster_scripts) repository for scripts and walkthrough for group members.

To setup a similar environment for your own group, install `decombinator` into an environment, download `Decombinator-Tags-FASTAs` to an accessible location, and adapt the scripts in the above repository to be compatible with your HPC service's job engine and directory structure.

Development Environment
-----------------------

To setup a development environment for `decombinator`, clone this repository to your computer:

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

How to collaborate
------------------

Please see our [code of conduct](./CODE_OF_CONDUCT.md) before contributing.

To report a bug or request support please post an issue [here](https://github.com/innate2adaptive/Decombinator/issues) outlining the problem faced and including error messages and logs where possible.

To propose a pull request, please create an issue first to discuss your proposed changes. We use [Black](https://github.com/psf/black) formatting, with a line length of 80.

For any support, commerical inquiries, or other requests, contact us [here](m.cowley@ucl.ac.uk).

Related reading
---------------

The published history of the development of the pipeline is covered in the following publications:

* [Thomas et al (2013), Bioinformatics: *Decombinator: a tool for fast, efficient gene assignment in T-cell receptor sequences using a finite state machine*](http://dx.doi.org/10.1093/bioinformatics/btt004)
* [Oakes et al (2017), Frontiers in Immunology: *Quantitative Characterization of the T Cell Receptor Repertoire of Naïve and Memory Subsets Using an Integrated Experimental and Computational Pipeline Which Is Robust, Economical, and Versatile*](https://doi.org/10.3389/fimmu.2017.01267)
* [Uddin et al (2019), Cancer Immunosurveillance: *An Economical, Quantitative, and Robust Protocol for High-Throughput T Cell Receptor Sequencing from Tumor or Blood*](http:/dx.doi.org/10.1007/978-1-4939-8885-3_2)
* [Uddin et al (2019), Methods in Enzymology: *Quantitative analysis of the T cell receptor repertoire*](https://doi.org/10.1016/bs.mie.2019.05.054)
* [Peacock et al (2020), Bioinformatics: *Decombinator V4: an improved AIRR compliant-software package for T-cell receptor sequence annotation*](https://doi.org/10.1093/bioinformatics/btaa758)
