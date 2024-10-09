##########################################################################
##################### WORKFLOW VARIABLES AND CODE ########################
##########################################################################

from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
from snakemake.utils import validate
from pathlib import Path
import sys
import os
import re
import pandas as pd
import polars as pl
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

sys.path.append((Path(workflow.basedir) / "../src/").as_posix())
sys.path.append((Path(workflow.basedir) / "scripts").as_posix())

nextpyper_version = "0.0.1"

from multi_seq_probes import group_probes, NoGrouping


report: "../report/workflow.rst"


SCHEMES_DIR = Path(workflow.basedir) / "schemes"

# Read inputs
probes_path = Path(config["args"]["probes"])
outdir = Path(config["args"]["output"])
path_samples = Path(config["args"]["input"])
pattern = config["args"]["probe_pattern"]
multi_probes = config["args"]["multi_probes"]
max_threads = config["args"]["threads"]

silva_db = Path(config["silva_db"])

## Read Workflow parameters:
pipeline = config["pipeline"]

# BBduk
bbduk_k = pipeline["matching_probes"]["bbduk_k"]
other_bbduk = pipeline["matching_probes"]["others"]

# MMseqs2
mmseq2_min_seq_id = pipeline["multi_probe_clustering"]["mmseq2_min_seq_id"]

# Spades
spades_k = "" if (argk := pipeline["spades"]["k"]) == "auto" else argk

# Split graph into probes
min_probe_cov = pipeline["split_graph_by_matching_probe"]["min_probe_coverage"]

# Blastx filtering
homolog_scf_min_cov = pipeline["blast_homolog_filtering"]["homolog_scf_min_cov"]
homolog_scf_min_idt = pipeline["blast_homolog_filtering"]["homolog_scf_min_idt"]

# Region separation
min_probe_contig_sim = pipeline["region_separation"]["min_probe_contig_sim"]
min_fragment_cov = pipeline["region_separation"]["min_fragment_cov"]
min_contig_length = pipeline["region_separation"]["min_contig_length"]

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


# Make useful structures for the inputs
sample_dict = SAMPLE_TABLE.set_index("sample").T.to_dict()
sample_list = list(sample_dict)

# Define pattern for matching final sequences
saute_seq_pattern = re.compile(
    r"(?P<sample>\w+)_Contig_(?P<probe>\w+)_(?P<cluster>\d+)-(?P<saute_info>[\d-]+)"
)


wildcard_constraints:
    sample=r"\w+",
    probe=r"\w+",
    cluster=r"\d+",
