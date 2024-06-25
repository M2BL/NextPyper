"""
Entrypoint for NextPiper

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

from pathlib import Path
import os
import click

from snaketool_utils.cli_utils import (
    OrderedCommands,
    run_snakemake,
    copy_config,
    echo_click,
)


def snake_base(rel_path):
    """Get the filepath to a Snaketool system file (relative to __main__.py)"""
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    """Read and print the version from the version file"""
    with open(snake_base("nextpiper.VERSION"), "r") as f:
        version = f.readline()
    return version


def print_citation():
    """Read and print the Citation information from the citation file"""
    with open(snake_base("nextpiper.CITATION"), "r") as f:
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
            default="nextpiper.out",
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
            default="nextpiper.log",
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
    nextpiper command --help"""
    pass


help_msg_extra = """
\b
CLUSTER EXECUTION:
nextpiper run ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           nextpiper run --input [file]
Specify threads:    nextpiper run ... --threads [threads]
Disable conda:      nextpiper run ... --no-use-conda 
Change defaults:    nextpiper run ... --snake-default="-k --nolock"
Add Snakemake args: nextpiper run ... --dry-run --keep-going --touch
Specify targets:    nextpiper run ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
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
    "--graph_simplification",
    "graph_simplification",
    help="Whether to simplying assembly graph after first assembly",
    is_flag=True,
    default=False,
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
    """Run NextPiper"""
    # Config to add or update in configfile
    merge_config = {"nextpiper": {"args": kwargs}}

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


cli.add_command(run)
cli.add_command(config)
cli.add_command(citation)


def main():
    cli()


if __name__ == "__main__":
    main()
