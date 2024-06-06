##########################################################################
##################### WORKFLOW VARIABLES AND CODE ########################
##########################################################################

from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
from snakemake.utils import validate
from pathlib import Path
import os
import pandas as pd

nextpiper_version = "0.0.1"


report: "../report/workflow.rst"


probes = Path(config["args"]["probes"])
outdir = Path(config["args"]["output"])
path_samples = config["args"]["_input"]
cols = ["sample_name", "path_forward", "path_reverse"]
SAMPLE_TABLE = pd.read_csv(path_samples, sep="\t", names=cols)
validate(SAMPLE_TABLE, schema="../schemes/sample_table.yaml")
sample_dict = SAMPLE_TABLE.set_index("sample_name").T.to_dict()
sample_list = list(sample_dict)
