from prefix_seqs import prefix_fasta, prefix_gfa
from gfa2fasta import paths2fasta

targets.append(outdir / "logs/dones/prefixing.done")


def get_graph_components(wildcards):
    w = wildcards
    checkpoint_output = checkpoints.split_graph_into_hmms.get(
        sample=wildcards.sample
    ).output[0]
    return Path(checkpoint_output) / f"{w.probe}.gfa"


rule prefix_components:
    input:
        get_graph_components,
    output:
        outdir / "assembled/prefixed/component_graphs/{sample}/{probe}.gfa",
    run:
        prefix = f"{wildcards.sample}-"
        prefix_gfa(input[0], output[0], prefix)


rule make_path_seqs:
    input:
        outdir / "assembled/prefixed/component_graphs/{sample}/{probe}.gfa",
    output:
        outdir / "assembled/prefixed/component_seqs/{sample}/{probe}.fasta",
    run:
        paths2fasta(Path(input[0]), Path(output[0]))


## See Rule 3.1 for further explanation
def aggregate_pref(wildcards):
    checkpoint_output = checkpoints.split_graph_into_hmms.get(
        sample=wildcards.sample
    ).output[0]
    global_match = glob_wildcards(Path(checkpoint_output) / "{probe}.gfa")
    # print(global_match)
    return expand(
        outdir
        / f"assembled/prefixed/component_seqs/{wildcards.sample}/{{probe}}.fasta",
        probe=global_match.probe,
    )


rule collect_prefixed_assemblies:
    input:
        aggregate_pref,
    output:
        chkpt=outdir / "logs/assembled/collect_pre/{sample}.chkpt",
    shell:
        "echo {input} | tr '[:space:]' '\n' >> {output.chkpt}"


# This checkpoint is necessary to reevaluate the DAG and make Snakemake
# realize that the prefixed asms exist and union_probes() can find them (4.1).
checkpoint done_prefix:
    input:
        chkpt=expand(
            outdir / "logs/assembled/collect_pre/{sample}.chkpt", sample=sample_list
        ),
    output:
        done=touch(outdir / "logs/dones/prefixing.done"),
