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
from Bio import SeqIO

sys.path.append((Path(workflow.basedir) / "../src/").as_posix())

nextpyper_version = "0.0.1"


report: "../report/workflow.rst"


# Parameters of the run
graph_simplification = config["args"]["graph_simplification"]

# Read inputs
probes = Path(config["args"]["probes"])
outdir = Path(config["args"]["output"])
path_samples = config["args"]["input"]

# Program configurations/parameters
# TAPER
taper_exec = Path(workflow.basedir) / config["taper_path"]
default_taper_params = config["taper_params_path"]

if not (path_taper_params := config["args"]["taper_params"]):
    path_taper_params = Path(workflow.basedir) / default_taper_params

# Trimal
trimal_gt = config["args"]["trimal_gt"]

# Validate Sample table
cols = ["sample_name", "path_forward", "path_reverse"]
SAMPLE_TABLE = pd.read_csv(path_samples, sep="\t", names=cols)
validate(SAMPLE_TABLE, schema="../schemes/sample_table.yaml")

# Validate probe names
probes_size = {probe.name: len(probe) for probe in SeqIO.parse(probes, "fasta")}
min_probe_size = min(probes_size.values())
probes_list = list(probes_size.keys())
PROBES = pd.DataFrame({"probe_name": probes_list})
validate(PROBES, schema="../schemes/probes.yaml")

# Make useful structures for the inputs
sample_dict = SAMPLE_TABLE.set_index("sample_name").T.to_dict()
sample_list = list(sample_dict)

# Define pattern for matching final sequences
saute_seq_pattern = re.compile(
    r"(?P<sample>\w+)_Contig_(?P<probe>\w+)_(?P<cluster>\d+)-(?P<saute_info>[\d-]+)"
)


wildcard_constraints:
    sample=r"\w+",
    probe=r"\w+",
    cluster=r"\d+",
