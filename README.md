# NextPiper
Recovery of homologous genes from targeted sequence capture data in higher ploidy samples.

---
# Dependencies
* [miniprot](https://bioconda.github.io/recipes/miniprot/README.html) ![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)
* [SPAdes](http://bioinf.spbau.ru/en/spades) 3.15.5
* [GetBlunted](https://github.com/vgteam/GetBlunted)
* [BubbleGun](https://github.com/fawaz-dabbaghieh/bubble_gun)
* [MAFFT](http://bioconda.github.io/recipes/mafft/README.html) ![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)
* [fastp](https://github.com/OpenGene/fastp)
* [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)


# Pipeline's structure
1. Read preprocessing: fastp
2. Read quality report (fastp, fastqc)
3. Target reads extraction: probe index building, read mapping (minimap/BWA2)
4. Contig assembly (SPAdes)
5. Prefixing assemblies (prefixing_seqs) Blunting graphs (get_blunted) and graph simplification (bubble_gun)
6. Homologs clustering (contig_cluster)
7. Alignment of clusters (MAFFT)
8. HMM profile building from alignments (hmm_build)
9. Graph exploration (graph_edge)
    1. with HMM profile (graph_edge)
    2. Minigraph (consensus needs to be extracted from HMM profile)
    3. SPAligner (works with assembly graphs, is it worth the time investment?)
10. Filtering of paths (to be implemented)
    optional phasing step
11. MSA per probe (MAFFT)
12. Annotate exon/intron regions
13. Summary statistics. 


# Output structure

```
Outdir/
├── QC_reads.json
├── Mapping_stats.tsv
├── Homolog_stats.tsv
└── Probe_alns
    ├── probe1_msa.fasta
    ├── probe2_msa.fasta
    ├── probe3_msa.fasta
    └── ...
```

# Installation

For now install the dependencies in a new conda environment. The installation of `nextpiper` will be contained in this environment. (**To Do:** improve deployment)

First clone the repository and install it.

```bash
git clone https://git.sorbus.ibot.cas.cz/m2b_ibot/nextpiper.git
cd nextpiper && pip install -e .
```

Now you can run `nextpiper` anywhere. 

## Running the test data

Now let us run a minimal test. Starting at the root directory of the repository:

```bash
# Make directories for the test and unpack the test data
mkdir -p minimal_test/test_data && cd minimal_test/test_data
tar xzf ../../nextpiper/data/test_data.tar.gz

# Prepare the sample table with the local paths
ls > ../inter.txt 
paste <(sed -E 's|(.*)_R[12].*|\1|' < ../inter.txt | uniq) <(grep "R1" ../inter.txt | xargs -I{} echo "$(pwd)/{}") <(grep "R2" ../inter.txt | xargs -I{} echo "$(pwd)/{}") > ../samples.tsv
cd ..

# Run nextpiper 
nextpiper run --input samples.tsv --probes ../../nextpiper/data/probes.fasta --output test_out -n
```

# Editing this README

When you're ready to make this README your own, just edit this file and use the handy template below (or feel free to structure it however you want - this is just a starting point!). Thanks to [makeareadme.com](https://www.makeareadme.com/) for this template.

## Suggestions for a good README

Every project is different, so consider which of these sections apply to yours. The sections used in the template are suggestions for most open source projects. Also keep in mind that while a README can be too long and detailed, too long is better than too short. If you think your README is too long, consider utilizing another form of documentation rather than cutting out information.

## Name
Choose a self-explaining name for your project.

## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.

## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection.

## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably include in the README.

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
