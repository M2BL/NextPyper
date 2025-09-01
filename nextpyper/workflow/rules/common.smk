##########################################################################
##################### WORKFLOW VARIABLES AND CODE ########################
##########################################################################

from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
from snakemake.utils import validate
from pathlib import Path
from collections import Counter, defaultdict
from operator import itemgetter
import json
import sys
import os
import re
import pandas as pd
import polars as pl
from more_itertools import last, one
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez

sys.path.append((Path(workflow.basedir) / "../src/").as_posix())

nextpyper_version = "0.0.1"

from multi_seq_probes import group_probes, NoGrouping


report: "../report/workflow.rst"


SCHEMES_DIR = Path(workflow.basedir) / "schemes"

# Read inputs
probes_path = Path(lookup("args/probes", within=config))
outdir = Path(lookup("args/output", within=config))
path_samples = Path(lookup("args/input", within=config))
probe_pattern = lookup("args/probe_pattern", within=config)
multi_probes = lookup("args/multi_probes", within=config)
max_threads = lookup("args/threads", within=config)
interseeds_use = lookup("args/interseeds", within=config)
reasm = lookup("args/reasm", within=config)

use_ref_cps = lookup("args/use_ref_cps", within=config)

blosum62 = Path(workflow.source_path(lookup("blosum62", within=config)))
cp_refs_map = Path(workflow.source_path(lookup("cp_refs_map", within=config)))

# Regex Patterns
SEED_PAT = lookup("regex_patterns/seed_pat", within=config)
SAUTE_PRE_FIX_PAT = lookup("regex_patterns/saute_pre_fix_pat", within=config)
SAUTE_POST_FIX_PAT = lookup("regex_patterns/saute_post_fix_pat", within=config)
TARGET_COLLAPSE_PAT = lookup("regex_patterns/saute_target_pat", within=config)

## Read Workflow parameters:
pipeline = lookup("pipeline", within=config)
scf_ext = lookup("scaffolds_extension", within=pipeline)
seeds_filt_params = lookup("homolog_filtering/seeds", within=pipeline)
homologs_filt_params = lookup("homolog_filtering/homologs", within=pipeline)
reg_sep = lookup("region_separation", within=pipeline)

# Validate Sample table
cols = ["sample", "path_forward", "path_reverse", "type"]
sample_table = pd.read_csv(path_samples, sep="\t", names=cols)
validate(sample_table, schema=(SCHEMES_DIR / "sample_table.yaml").resolve())
sample_list = sample_table["sample"].to_list()

# Validate probes
probes = SeqIO.to_dict(SeqIO.parse(probes_path, "fasta"))
probes_list = list(probes.keys())
PROBES = pd.DataFrame({"probe_name": probes_list})
validate(PROBES, schema=(SCHEMES_DIR / "probes.yaml").resolve())

# Check grouping of the probe set
# Multi-seq probe set
if multi_probes:
    try:
        probe_hier = {
            probe: [rec.id for rec in recs]
            for probe, recs in group_probes(
                list(probes.values()), probe_pattern
            ).items()
        }
        probes_list = list(probe_hier.keys())
    except NoGrouping as e:
        e.add_note(
            "Either this is a single-probe set or the pattern is not appropiate for multi-probe"
        )
        raise
else:
    pat = re.compile(probe_pattern)
    try:
        probes_list = [pat.search(probe)[1] for probe in probes]
    except TypeError:
        names = [probe.id for probe in probes if not pat.search(probe.id)]
        raise Exception(
            f"Error: pattern {probe_pattern} does not match all probes:\n{"\n".join(names)}"
        )


wildcard_constraints:
    sample=r"\w+",
    probe=r"\w+",
    cluster=r"\d+",
