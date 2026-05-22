"""
Entrypoint for NextPyper
"""

from pathlib import Path
import sys
import tempfile

import yaml
import rich_click as click
from snaketool_utils.cli_utils import run_snakemake, copy_config, echo_click

sys.path.append(str((Path(__file__).parent / "workflow/scripts").resolve()))
sys.path.append(str((Path(__file__).parent / "src").resolve()))

from sample_table import make_table
from multi_seq_probes import check_probes
from summarize_results import summarize_workflow


def snake_base(rel_path):
    """Get the filepath to a Snaketool system file (relative to __main__.py)"""
    return (Path(__file__).parent / rel_path).resolve()


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
        return str(Path(ctx.params["output"]) / value)
    return value


# Load CLI configuration.
config_path = Path(__file__).parent / "config/rich_config.yaml"
if config_path.exists():
    # click.rich_click.__dict__.update(yaml.safe_load(config_path.read_text()))
    config_data = yaml.safe_load(config_path.read_text())
    for key, value in config_data.items():
        # Version 1.9.x uses uppercase for global variables
        setattr(click.rich_click, key.upper(), value)


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
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
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default=str(snake_base(Path("workflow") / "conda")),
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
            default=snake_base(Path("config/config.yaml")),
            hidden=True,
        ),
        click.argument("snake_args", nargs=-1),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(
    cls=click.RichGroup, context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option(get_version(), "-v", "--version", is_flag=True)
def cli():
    """Recovery of homoeologous loci from target capture data in higher ploidy samples.
    \n
    For more options, run:
    nextpyper command --help"""
    pass


# ToDo: Refine this help message
gather_msg = """
The data directory is expected to have raw paired reads
files per each sample (forward, reverse).
"""

# ToDo: Refine this help message
validate_msg = """

"""

# ToDo: Refine this help message
summarize_msg = """

"""

prepare_msg = """
By default, the environments will be created in NextPyper's installation folder. 
Use --conda-prefix to specify a custom location.
"""


@click.command(
    # epilog=help_msg_extra,
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
    "--output",
    help="Output directory",
    type=click.Path(dir_okay=True, writable=True, readable=True),
    default="nextpyper.out",
    show_default=True,
)
@click.option(
    "--probes",
    "probes",
    help="Probes used in the experiment (fasta)",
    type=click.Path(readable=True, exists=True),
    required=True,
)
@click.option(
    "--pattern",
    "probe_pattern",
    help="For multi-probes, pattern used to group the probe sequences (RegEx). It needs one capture group.",
    type=str,
    default=r"(.*)",
    show_default=True,
)
@click.option(
    "--multi-probes/--single-probes",
    "multi_probes",
    help="Whether the probe set has multiple or a single sequence per probe ",
    default=True,
    show_default=True,
)
@click.option(
    "--use-ploidy/--no-ploidy",
    "use_ploidy",
    help="""Whether to use the ploidy information from the samples. An extra
            ploidy column is expected in the input sample table. Ploidy is 
            used to inform the expected number of homologs.""",
    default=False,
    show_default=True,
)
@click.option(
    "--interseeds",
    "interseeds",
    help="""Which inter-sample seeds to use for SAUTE assembly. Interseeds 
            help to boost probe recovery at the expense of higher 
            computation time during assembly. The default 'sister' infers 
            sister-samples that are likely the most informative to 
            get inter-sample seeds from.""",
    type=click.Choice(("all", "sister", "none")),
    default="sister",
    show_default=True,
)
@click.option(
    "--reasm-complex-probes/--no-reasm",
    "reasm",
    help="""Whether to reassemble the most complex ("explosive") probes found in
            each sample. This second assembly is tailored to better resolve such 
            complexity.""",
    default=True,
    show_default=True,
)
@click.option(
    "--use-ref-cps/--no-ref-cps",
    "use_ref_cps",
    help="""Download reference chloroplasts for cpDNA filtering.""",
    default=True,
    show_default=True,
)
@click.option(
    "--ref-cps",
    "ref_cps",
    help="""Comma separated list of reference chloroplasts to download for
            cp filtering. Check cps_seqids.csv for a full list of values.""",
    type=str,
    default="Ambtr,Arath,Orysa,UMUL,RQNK",
    show_default=True,
)
@click.option(
    "--custom-cps",
    "custom_cps",
    help="""Custom cps to use for cpDNA filtering (fasta). 
            These cps are used in addition to the ref cps.""",
    type=click.Path(readable=True, exists=True),
    required=False,
    default=None,
)
@common_options
def run(**kwargs):
    """Run NextPyper"""
    # Config to add or update in configfile
    merge_config = {"nextpyper": {"args": kwargs}}

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(Path("workflow/Snakefile")),
        merge_config=merge_config,
        **kwargs
    )


@click.command()
@click.option(
    "--output",
    help="Output directory",
    type=click.Path(dir_okay=True, writable=True, readable=True),
    default="nextpyper.out",
    show_default=True,
)
# @common_options
def config(configfile, system_config, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile, system_config=system_config)


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


@click.option(
    "--output",
    help="Output sample table",
    type=click.Path(writable=True, readable=True),
    default="sample.tsv",
    show_default=True,
)
@click.option(
    "--input",
    "input",
    help="Path to data directory",
    type=click.Path(readable=True, exists=True),
    required=True,
)
@click.command(epilog=gather_msg)
def gather(**kwargs):
    """Generate a sample table given a data directory"""

    datadir = Path(kwargs["input"])
    outfile = Path(kwargs["output"])

    make_table(datadir, outfile)


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
@click.option(
    "--pattern",
    help="Pattern used to group the probes (RegEx)",
    type=str,
    default=r"(.*)",
    show_default=True,
)
@click.option(
    "--probes",
    "probes",
    help="Path to probes files",
    type=click.Path(readable=True, exists=True),
    required=True,
)
@click.command(epilog=validate_msg)
def validate(**kwargs):
    """Validate a probes file for running NextPyper"""

    probes_path = Path(kwargs["probes"])
    pattern = kwargs["pattern"]
    outsummary = Path(kwargs["output"]) if kwargs["output"] else None
    outhierarchy = Path(kwargs["hierarchy"]) if kwargs["hierarchy"] else None

    check_probes(probes_path, pattern, outsummary, outhierarchy)


@click.option(
    "--output",
    "output",
    help="Output summary table",
    type=click.Path(writable=True, path_type=Path),
    default="run_stats.csv",
    show_default=True,
)
@click.option(
    "--rundir",
    "rundir",
    help="Path to run directory",
    type=click.Path(readable=True, exists=True, path_type=Path),
    required=True,
)
@click.command(epilog=summarize_msg)
def summarize(**kwargs):
    """Summarize the results of a NextPyper run"""

    run_directory_path = kwargs["rundir"]
    out_table_path = kwargs["output"]
    tab_file = kwargs["seqs_per_probe"]

    # df, table = summarize_workflow(run_directory_path)
    df = summarize_workflow(run_directory_path)
    df.to_csv(out_table_path, index=False)

    # if tab_file:
    #     table.T.to_csv(tab_file, float_format="%.2f")


@click.option(
    "--conda-prefix",
    default=str(snake_base(Path("workflow") / "conda")),
    help="Path where to put the conda environments",
    type=click.Path(),
    show_default=False,
)
@click.command(epilog=prepare_msg)
def prepare(**kwargs):
    """Prepare the conda environments to run NextPyper"""

    mock_dir = Path(__file__).parent / "data/mock"
    kwargs.update(yaml.safe_load((mock_dir / "mock_args.yaml").read_text()))

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        kwargs["system_config"] = str(Path(__file__).parent / "config/config.yaml")
        kwargs["input"] = str(mock_dir / "sample.tsv")
        kwargs["probes"] = str(mock_dir / "probes.fasta")
        kwargs["configfile"] = str(tmpdir / "config.yaml")
        kwargs["log"] = str(tmpdir / "nextpyper.log")
        kwargs["output"] = str(tmpdir)

        # Config to add or update in configfile
        merge_config = {"nextpyper": {"args": kwargs}}

        run_snakemake(
            # Full path to Snakefile
            snakefile_path=snake_base(Path("workflow/Snakefile")),
            merge_config=merge_config,
            **kwargs
        )


cli.add_command(run)
cli.add_command(gather)
cli.add_command(validate)
cli.add_command(summarize)
cli.add_command(prepare)
cli.add_command(config)
cli.add_command(citation)


def main():
    cli()


if __name__ == "__main__":
    main()
