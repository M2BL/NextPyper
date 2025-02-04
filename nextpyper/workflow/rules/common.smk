##########################################################################
##################### WORKFLOW VARIABLES AND CODE ########################
##########################################################################

from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
from snakemake.utils import validate
from pathlib import Path
from collections import Counter, defaultdict
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

# MMseqs2
mmseq2_min_seq_id = lookup("multi_probe_clustering/mmseq2_min_seq_id", within=pipeline)

# Spades
spades_k = "" if (argk := lookup("spades/k", within=pipeline)) == "auto" else argk

# MMseqs prefiltering
mmseq_prefilt_sens = lookup("mmseqs_prefiltering/sensitivity", within=pipeline)

# Split graph into probes
min_probe_cov = lookup(
    "split_graph_by_matching_probe/min_probe_coverage", within=pipeline
)

# Scaffold extension
floor_len_extension = lookup(
    "scaffolds_extension/exploration/floor_len_extension", within=pipeline
)
plen_scaling_factor = lookup(
    "scaffolds_extension/exploration/probe_len_scaling", within=pipeline
)

# MMseqs matching
mmseq_fields = lookup("mmseqs_matching/fields", within=pipeline)
mmseq_evalue = lookup("mmseqs_matching/evalue", within=pipeline)
min_orf_len = lookup("mmseqs_matching/min_orf_len", within=pipeline)
mmseq_sens = lookup("mmseqs_matching/sensitivity", within=pipeline)

# MMseqs filtering
homolog_scf_min_cov = lookup("homolog_filtering/homolog_scf_min_cov", within=pipeline)
homolog_scf_min_idt = lookup("homolog_filtering/homolog_scf_min_idt", within=pipeline)

# Region separation
min_probe_contig_sim = lookup("region_separation/min_probe_contig_sim", within=pipeline)
min_fragment_cov = lookup("region_separation/min_fragment_cov", within=pipeline)
min_contig_length = lookup("region_separation/min_contig_length", within=pipeline)

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
    probes_list = [cpat.search(probe)[1] for probe in all_probes_list]


# Make useful structures for the inputs
sample_dict = SAMPLE_TABLE.set_index("sample").T.to_dict()
sample_list = list(sample_dict)


wildcard_constraints:
    sample=r"\w+",
    probe=r"\w+",
    cluster=r"\d+",
