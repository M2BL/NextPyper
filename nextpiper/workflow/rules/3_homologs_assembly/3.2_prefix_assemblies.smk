from prefix_seqs import prefix_fasta, prefix_gfa

targets.append(outdir / "logs/dones/prefixing.done")


rule prefix_assemblies:
    input:
        contigs=outdir / "assembled/spades/{sample}/{probe}/contigs.fasta",
        gfa=outdir
        / "assembled/spades/{sample}/{probe}/assembly_graph_with_scaffolds.gfa",
    output:
        contigs=outdir / "assembled/prefixed/{sample}/{probe}/contigs.fasta",
        gfa=outdir
        / "assembled/prefixed/{sample}/{probe}/assembly_graph_with_scaffolds.gfa",
    log:
        outdir / "logs/assembled/spades/prefixing/{sample}/{probe}.log",
    threads: 1
    run:
        prefix = f"{wildcards.sample}-{wildcards.probe}-"
        prefix_fasta(input.contigs, output.contigs, prefix)
        prefix_gfa(input.gfa, output.gfa, prefix)


## See Rule 3.1 for further explanation
def aggregate_pref(wildcards):
    checkpoint_output = checkpoints.distribute_reads.get(**wildcards).output[0]
    return expand(
        outdir / "assembled/prefixed/{sample}/{probe}/contigs.fasta",
        sample=wildcards.sample,
        probe=glob_wildcards(os.path.join(checkpoint_output, "{probe}.bam")).probe,
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
