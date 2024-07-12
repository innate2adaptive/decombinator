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
TCR repertoire sequencing (TCRseq) offers a powerful means to investigate biological samples to observe sequence distribtuions.
However current high-throughput sequencing (HTS) technologies can produce large amounts of raw data which will also unavoidably contain errors relative to the original input molecules.

Decombinator addresses the problem of large datasets through speed - employing a rapid and [highly efficient string matching algorithm](https://figshare.com/articles/Aho_Corasick_String_Matching_Video/771968) to search the FASTQ files produced by HTS machines for rearranged TCR sequences.
The central algorithm searches for 'tag' sequences, the presence of which uniquely indicates the inclusion of particular V or J genes in a recombination.
If V and J tags are found, Decombinator can then deduce where the ends of the germline V and J gene sections are (i.e. how much nucleotide removal occurred during V(D)J recombination), and what nucleotide sequence (the 'insert sequence') remains between the two.
These five pieces of information - the V and J genes used, how many deletions each had and the insert sequence - contain all of the information required to reconstruct the whole TCR nucleotide sequence, in a more readily stored and analysed way.
The pipeline also handles FASTQ data with many reads of the same molecule due to over-amplification.
It does this by using unique molecular identifiers to "collapse" down duplicate reads.

Running decombinator is easy. All that is required is HTS read data and a few arguments specifying the chemistry of your molecules:

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

In this mode, all three main components of the pipeline are applied to the data sequentially: `decombine`, `collapse`, and `translate`.
For full detail of each of these components please read the [documentation](tbd).

If instead, you wish to run just one component of the pipeline, these are all directly accessible via their respective sub-parser:

```shell
decombinator decombine -in XXXX.fq -c b -br R2 -bl 42
```

```shell
decombinator collapse -in XXXX.n12 -c b -ol M13
```

```shell
decombinator translate -in XXXX.freq -c b -ol M13
```

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
pip install .
```

Tests can be then run via:

```shell
pytest
```

How to collaborate
------------------

To report a bug please post an issue [here](https://github.com/innate2adaptive/Decombinator/issues).

Related reading
---------------

The published history of the development of the pipeline is covered in the following publications:

* [Thomas et al (2013), Bioinformatics: *Decombinator: a tool for fast, efficient gene assignment in T-cell receptor sequences using a finite state machine*](http://dx.doi.org/10.1093/bioinformatics/btt004)
* [Oakes et al (2017), Frontiers in Immunology: *Quantitative Characterization of the T Cell Receptor Repertoire of Naïve and Memory Subsets Using an Integrated Experimental and Computational Pipeline Which Is Robust, Economical, and Versatile*](https://doi.org/10.3389/fimmu.2017.01267)
* [Uddin et al (2019), Cancer Immunosurveillance: *An Economical, Quantitative, and Robust Protocol for High-Throughput T Cell Receptor Sequencing from Tumor or Blood*](http:/dx.doi.org/10.1007/978-1-4939-8885-3_2)
* [Uddin et al (2019), Methods in Enzymology: *Quantitative analysis of the T cell receptor repertoire*](https://doi.org/10.1016/bs.mie.2019.05.054)
* [Peacock et al (2020), Bioinformatics: *Decombinator V4: an improved AIRR compliant-software package for T-cell receptor sequence annotation*](https://doi.org/10.1093/bioinformatics/btaa758)

