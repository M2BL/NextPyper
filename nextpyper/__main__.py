"""
Entrypoint for NextPyper

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

from pathlib import Path
import os
import click
import sys

from snaketool_utils.cli_utils import (
    OrderedCommands,
    run_snakemake,
    copy_config,
    echo_click,
)

sys.path.append(str((Path(__file__).parent / "workflow/scripts").resolve()))
from sample_table import make_table
from multi_seq_probes import check_probes


def snake_base(rel_path):
    """Get the filepath to a Snaketool system file (relative to __main__.py)"""
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    """Read and print the version from the version file"""
    with open(snake_base("nextpyper.VERSION"), "r") as f:
        version = f.readline()
    return version


def print_citation():
    """Read and print the Citation information from the citation file"""
    with open(snake_base("nextpyper.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def default_to_output(ctx, param, value):
    """Callback for click options; places value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "--output",
            help="Output directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="nextpyper.out",
            show_default=True,
        ),
        click.option(
            "--configfile",
            default="config.yaml",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/config.yaml]",
        ),
        click.option(
            "--threads", help="Number of threads to use", default=1, show_default=True
        ),
        click.option(
            "--profile",
            default=None,
            help="Snakemake profile to use",
            show_default=False,
        ),
        click.option(
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default=snake_base(os.path.join("workflow", "conda")),
            help="Custom conda env directory",
            type=click.Path(),
            show_default=False,
        ),
        click.option(
            "--snake-default",
            multiple=True,
            default=[
                "--printshellcmds",
                "--nolock",
                "--show-failed-logs",
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        ),
        click.option(
            "--log",
            default="nextpyper.log",
            callback=default_to_output,
            hidden=True,
        ),
        click.option(
            "--system-config",
            default=snake_base(os.path.join("config", "config.yaml")),
            hidden=True,
        ),
        click.argument("snake_args", nargs=-1),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(
    cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option(get_version(), "-v", "--version", is_flag=True)
def cli():
    """Recovery of homologous genes from targeted sequence capture data for higher ploidy samples
    \b
    For more options, run:
    nextpyper command --help"""
    pass


help_msg_extra = """
\b
CLUSTER EXECUTION:
nextpyper run ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           nextpyper run --input [file]
Specify threads:    nextpyper run ... --threads [threads]
Disable conda:      nextpyper run ... --no-use-conda 
Change defaults:    nextpyper run ... --snake-default="-k --nolock"
Add Snakemake args: nextpyper run ... --dry-run --keep-going --touch
Specify targets:    nextpyper run ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""

# ToDo: Refine this help message
sample_table_msg = """
The data directory is expected to have raw paired reads
files per each sample (forward, reverse).
"""

# ToDo: Refine this help message
validate_probes_msg = """

"""


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
    "--input",
    "input",
    help="Input sample table",
    type=click.Path(readable=True, exists=True),
    required=True,
)
@click.option(
    "--probes",
    "probes",
    help="Probes used in the experiment (fasta)",
    type=click.Path(readable=True, exists=True),
    required=True,
)
@click.option(
    "--multi-probes/--single-probes",
    "multi_probes",
    help="Whether the probe set has multiple or a single sequence per probe ",
    default=True,
    show_default=True,
)
@click.option(
    "--taper_parameters",
    "taper_params",
    help="Parameters file to use when running TAPER (-p). See TAPER's docs.",
    type=click.Path(readable=True, exists=True),
)
@click.option(
    "--trimal_gt",
    "trimal_gt",
    help="1 - (fraction of sequences with a gap allowed)in Trimal.",
    type=float,
    default=0.2,
)
@common_options
def run(**kwargs):
    """Run NextPyper"""
    # Config to add or update in configfile
    merge_config = {"nextpyper": {"args": kwargs}}

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "Snakefile")),
        merge_config=merge_config,
        **kwargs
    )


@click.command()
@common_options
def config(configfile, system_config, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile, system_config=system_config)


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


@click.option(
    "--input",
    "input",
    help="Path to data directory",
    type=click.Path(readable=True, exists=True),
    required=True,
)
@click.option(
    "--output",
    help="Output sample table",
    type=click.Path(writable=True, readable=True),
    default="sample.tsv",
    show_default=True,
)
@click.command(epilog=sample_table_msg)
def make_sample_table(**kwargs):
    """Generate a sample table given a data directory"""

    datadir = Path(kwargs["input"])
    outfile = Path(kwargs["output"])

    make_table(datadir, outfile)


@click.option(
    "--probes",
    "probes",
    help="Path to probes files",
    type=click.Path(readable=True, exists=True),
    required=True,
)
@click.option(
    "--pattern",
    help="Pattern used to group the probes (RegEx)",
    type=str,
    default=r"(\d{4})$",
    # click.Path(writable=True, readable=True),
    # default="sample.tsv",
    show_default=True,
)
@click.option(
    "--write_hierarchy",
    "hierarchy",
    help="Write grouping hierarchy to this file",
    type=click.Path(writable=True, readable=True),
    default="",
    # show_default=True,
)
@click.option(
    "--write_summary",
    "output",
    help="Write summary of grouping to this file",
    type=click.Path(writable=True, readable=True),
    default="",
    # show_default=True,
)
@click.command(epilog=validate_probes_msg)
def validate_probes(**kwargs):
    """Validate a probes file for running NextPyper"""

    probes_path = Path(kwargs["probes"])
    pattern = kwargs["pattern"]
    outsummary = Path(kwargs["output"]) if kwargs["output"] else None
    outhierarchy = Path(kwargs["hierarchy"]) if kwargs["hierarchy"] else None

    check_probes(probes_path, pattern, outsummary, outhierarchy)


cli.add_command(run)
cli.add_command(make_sample_table)
cli.add_command(validate_probes)
cli.add_command(config)
cli.add_command(citation)


def main():
    cli()


if __name__ == "__main__":
    main()
