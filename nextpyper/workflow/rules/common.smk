##########################################################################
##################### WORKFLOW VARIABLES AND CODE ########################
##########################################################################

from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
from snakemake.utils import validate
from pathlib import Path
from collections import Counter, defaultdict
import json
import sys
import os
import re
import pandas as pd
import polars as pl
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez

sys.path.append((Path(workflow.basedir) / "../src/").as_posix())
sys.path.append((Path(workflow.basedir) / "scripts").as_posix())

nextpyper_version = "0.0.1"

from prefix_seqs import prefix_fasta
from multi_seq_probes import group_probes, NoGrouping


report: "../report/workflow.rst"


SCHEMES_DIR = Path(workflow.basedir) / "schemes"

# Read inputs
probes_path = Path(lookup("args/probes", within=config))
outdir = Path(lookup("args/output", within=config))
path_samples = Path(lookup("args/input", within=config))
pattern = lookup("args/probe_pattern", within=config)
multi_probes = lookup("args/multi_probes", within=config)
max_threads = lookup("args/threads", within=config)

blosum62 = workflow.source_path(config["blosum62"])

## Read Workflow parameters:
pipeline = lookup("pipeline", within=config)

# Fastp
trim_qual = lookup("preprocessing/trim_qual", within=pipeline)
trim_min_len = lookup("preprocessing/min_len", within=pipeline)

# Spades
spades_k = lookup("spades/k", within=pipeline)

# Scaffold extension
scf_ext = lookup("scaffolds_extension", within=pipeline)
floor_len_extension = lookup("exploration/floor_len_extension", within=scf_ext)
ceil_len_extension = lookup("exploration/ceiling_len_extension", within=scf_ext)
plen_scaling_factor = lookup("exploration/probe_len_scaling", within=scf_ext)
max_intron_size = lookup("exploration/max_intron_size", within=scf_ext)
max_extensions = lookup("output/max_extensions", within=scf_ext)

# MMseqs matching
mmseq_fields = lookup("mmseqs_matching/fields", within=pipeline)
mmseq_evalue = lookup("mmseqs_matching/evalue", within=pipeline)
min_orf_len = lookup("mmseqs_matching/min_orf_len", within=pipeline)
mmseq_sens = lookup("mmseqs_matching/sensitivity", within=pipeline)

# MMseqs filtering
homolog_scf_min_cov = lookup("homolog_filtering/homolog_scf_min_cov", within=pipeline)
homolog_scf_min_idt = lookup("homolog_filtering/homolog_scf_min_idt", within=pipeline)

# Region separation
reg_sep = lookup("region_separation", within=pipeline)
min_probe_scaffold_sim = lookup("min_probe_scaffold_sim", within=reg_sep)
min_fragment_cov = lookup("min_fragment_cov", within=reg_sep)
min_exonic_length = lookup("min_exonic_length", within=reg_sep)

# Validate Sample table
cols = ["sample", "path_forward", "path_reverse", "type"]
SAMPLE_TABLE = pd.read_csv(path_samples, sep="\t", names=cols)
validate(SAMPLE_TABLE, schema=(SCHEMES_DIR / "sample_table.yaml").resolve())

# Validate probes
probes = list(SeqIO.parse(probes_path, "fasta"))
probes_size = {probe.id: len(probe) for probe in probes}
min_probe_size = min(list(probes_size.values()))
probes_list = list(probes_size.keys())
PROBES = pd.DataFrame({"probe_name": probes_list})
validate(PROBES, schema=(SCHEMES_DIR / "probes.yaml").resolve())

# Check grouping of the probe set
# Multi-seq probe set
if multi_probes:
    try:
        probe_hier = {
            probe: [rec.id for rec in recs]
            for probe, recs in group_probes(probes, pattern).items()
        }
        all_probes_list = probes_list
        probes_list = list(probe_hier.keys())
    except NoGrouping:
        # Single-sequence probe set
        multi_probes = False
else:
    cpat = re.compile(pattern)
    all_probes_list = probes_list
    try:
        probes_list = [cpat.search(probe)[1] for probe in all_probes_list]
    except TypeError:
        names = [probe.id for probe in probes if not cpat.search(probe.id)]
        raise Exception(
            f"Error: pattern {pattern} does not match all probes:\n{"\n".join(names)}"
        )


# Make useful structures for the inputs
sample_dict = SAMPLE_TABLE.set_index("sample").T.to_dict()
sample_list = list(sample_dict)


wildcard_constraints:
    sample=r"\w+",
    probe=r"\w+",
    cluster=r"\d+",
