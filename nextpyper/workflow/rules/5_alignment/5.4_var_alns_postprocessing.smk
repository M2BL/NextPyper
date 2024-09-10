use rule taper as taper_vars with:
    input:
        outdir / "var_aligned/var_alns/{probe}.fasta",
    output:
        outdir / "var_aligned/taper/{probe}.fasta",
    log:
        outdir / "logs/var_aligned/taper/{probe}.log",


use rule trimal as trimal_vars with:
    input:
        outdir / "var_aligned/taper/{probe}.fasta",
    output:
        outdir / "var_aligned/trimal/{probe}.fasta",
    log:
        outdir / "logs/var_aligned/trimal/{probe}.log",


## Is this collections necessary?
# def aggregate_clean_var_alns(wildcards):
#     checkpoint_output = checkpoints.group_variants_by_probe.get(
#         probe=wildcards.probe
#     ).output[0]
#     glob_match = glob_wildcards(Path(checkpoint_output) / f"{wildcards.probe}.fasta")

#     return expand(
#         outdir / "var_aligned/trimal/{probe}.fasta",
#         probe=wildcards.probe,
#     )


rule collect_clean_var_alignments:
    input:
        # aggregate_clean_var_alns,
        outdir / "var_aligned/trimal/{probe}.fasta",
    output:
        chkpt=outdir / "logs/var_aligned/collect_clean_var_alns/{probe}.chkpt",
    shell:
        "echo {input} | tr '[:space:]' '\n' >> {output.chkpt}"


checkpoint done_var_cleaning:
    input:
        chkpt=expand(
            outdir / "logs/var_aligned/collect_clean_var_alns/{probe}.chkpt",
            probe=probes_list,
        ),
    output:
        done=touch(outdir / "logs/dones/clean_var_alns.done"),
