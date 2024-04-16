##########################################################################
##################### WORKFLOW VARIABLES AND CODE ########################
##########################################################################

from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
from snakemake.utils import validate
import os

nextpiper_version = "0.0.1"

configfile: "config/config.yaml"

report: "../report/workflow.rst"
