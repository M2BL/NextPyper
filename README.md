# NextPyper
Recovery of homoeologous loci from targeted capture data in higher ploidy samples.

## Table of contents

- [NextPyper](#nextpyper)
  - [Table of contents](#table-of-contents)
  - [Installation](#installation)
    - [Conda/Mamba](#condamamba)
    - [Pip](#pip)
    - [Manual](#manual)
  - [Usage](#usage)
    - [1. Preparing your probes](#1-preparing-your-probes)
      - [Single-probe sets](#single-probe-sets)
      - [Multi-probe sets](#multi-probe-sets)
    - [2. Preparing your samples](#2-preparing-your-samples)
    - [3. Running NextPyper](#3-running-nextpyper)
  - [Output](#output)


## Installation

Regardless of how you install it, `NextPyper` requires `conda` to be available in your PATH, to handle and deploy the environments to run the workflow.

### Conda/Mamba

If you already have [bioconda](https://bioconda.github.io/) set up, you can create a new environment and simply run:

```bash
conda install nextpyper # or mamba
```

### Pip

Alternatively you can use pip:

```
pip install nextpyper
```

### Manual 

You can also install NextPyper manually with pip. We recommend to create a new environment.

```bash
# from you new environment with pip
git clone git@github.com:M2BL/NextPyper.git
cd NextPyper
pip install .
```

## Usage 

### 1. Preparing your probes

NextPyper expects a nucleotide fasta file with the sequences of the targeted loci, what we call the probe set. Probe sets can be divided in two categories:

- **Multi probe sets**: Contains multiple sequences per locus/gene (*e.g.* Angiosperms353).
- **Single probe sets**: Contains a single sequence per locus. Usually, custom probe sets designed from a single reference genome are in this category.

Sequence names in the probe set are required to be simple (letters, numbers and underscore, characters only). Additionally no duplicated names are allowed.

Currently, only coding regions are supported.

#### Single-probe sets

If you are using a single-probe set, only need to worry sticking to the naming conventions described above. To check your probes, you can use `nextpyper validate`:

```bash
$ nextpyper validate --probes probes.fasta 

# Probe sequence names comply with naming convention.
# Pattern (.*) yielded no grouping of the sequences.
# Either this is a single-probe set or the pattern is not appropiate for multi-probe mode
```

If you get the message above, your probe set is well formatted to run in single-probes mode.

#### Multi-probe sets

If using a multi-probe set, NextPyper needs to know which sequences target the same locus. The option `--pattern` expects a RegEx (with one capture group) that will inform about that hierarchy based on the name of the probes. 

For instance, for the Angiosperms353 probe set with sequence names:

```
AJFN_probe4471
Ambtr_probe4471
BERS_probe4527
ZENX_probe4527
Arath_probe4691
QUTB_probe4691
```

We can use `--pattern "(\d+)$"`, which will group the sequences by the number at the end of the sequence name (which is the LocusID).

You can use the `validate` subcommand to check if your pattern groups your probes:

```bash
$ nextpyper validate --probes kew_probes.fasta --pattern "(\d+)$"

# Probe sequence names comply with naming convention.
# The pattern yields 353 groups from 4781 probe sequences
```

Use the options `--out_summary` and `out_hierarchy` to get more details about the grouping.

### 2. Preparing your samples

NextPyper expects a 4-column table that defines the identity and ploidy of your samples, and where to find their associated data.

The table requires a header with 4 columns (see an example here):

- sample: the name of the sample.
- path_forward: absolute path to the forward reads. It can be .gz compressed.
- path_reverse: absolute path to the reverse reads. It can be .gz compressed.
- ploidy: the ploidy of the sample. If 0, ploidy won't be taken into account to process the sample.

You can use the `gather` subcommand to help you create the sample table. From a data directory ([like this one](https://github.com/M2BL/NextPyper/tree/main/nextpyper/data/mock/reads)) with only the paired-end data (forward and reverse files) you can run:

```bash
nextpyper gather --input data_folder --output sample.tsv
```

To get the input the sample table for your samples. See an example of the result [here](https://github.com/M2BL/NextPyper/blob/main/nextpyper/data/mock/sample.tsv).

The sample names are inferred based on the common prefix between the forward and reverse files. 

### 3. Running NextPyper

Finally, you can simply run NextPyper, for single-probes:

```bash
nextpyper run --threads 16 --single-probes --input input_samples.tsv --probes custom_probes.fasta --output out_single
```

Or for multi-probes, specifying the appropiate pattern:

```bash
nextpyper run --threads 16 --input new_hier_samples.tsv --probes kew_probes.fasta --pattern "(\d+)$" --output out_multi
```


## Output

The main output of NextPyper are the samples' homologous sequences for each tageted loci. We generate 3 versions of these sequences, depending on the final postprocesing:

- Exons: includes only the exonic part of the sequence.
- Genetigs: includes exons and introns.
- Supercontigs: The full sequence, including exons, introns and flanking regions. 

You may want to explore which type of sequence suits better your purposes.

NextPyper produces several outputs you might be interested in:

- Homologous sequences per sample: `<outdir>/homolog_prospection/region_separation/consolidated/per_sample`
- Homologous sequences per loci: `<outdir>/homolog_prospection/region_separation/separation_output/scfs`
- Alignments: `<outdir>/homolog_prospection/region_separation/alns`
