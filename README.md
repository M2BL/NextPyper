# NextPyper

**Recovery of homoeologous loci from target capture data in higher ploidy samples.**

`NextPyper` is a bioinformatics pipeline specifically designed to disentangle and recover distinct homoeologs, paralogs, and even ohnologs from polyploid target capture datasets.

## Key Features

- **Built for polyploids**: Handles polyploid samples across a wide range of ploidal levels (from diploids to 10x+ polyploids) without requiring high-quality parental reference genomes.
- **Uncompromised sensitivity**: Strives to recover all distinct homologs (homoeologs, paralogs, and even ohnologs) present in your data.
- **Inter-sample seed borrowing**: Rescues targeted loci in low-coverage or diverged samples by leveraging homologous sequences from related samples.
- **Automated graph complexity resolution**: Automatically identifies dense assembly graph components ("tribbles"), re-assembles complex loci with adaptive parameters, and flags putative chimeric sequences.
- **Streamlined workflow**: Simple command-line interface powered by [Snakemake](https://snakemake.github.io/), offering automatic environment management, parallel processing, and seamless execution.

## Documentation

Please refer to the [official documentation](https://m2bl.github.io/NextPyper/), which is more complete.

## Installation

`NextPyper` is available in [PyPI](https://pypi.org/project/NextPyper/), and soon in bioconda.

Regardless of the installation method, `NextPyper` requires `conda` in your PATH, to manage and deploy workflow environments.

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
# from your new environment with pip
git clone git@github.com:M2BL/NextPyper.git
cd NextPyper
pip install .
```

## License
MIT