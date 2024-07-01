targets.append(outdir / "logs/dones/clean_alns.done")


rule taper:
    input:
        outdir / "aligned/cluster_alns/{probe}/{probe}_{cluster}.fasta",
    output:
        outdir / "aligned/taper/{probe}/{probe}_{cluster}.fasta",
    params:
        str(taper_exec) + " -m - -a - -p " + str(path_taper_params),
    log:
        outdir / "logs/aligned/taper/{probe}/{probe}_{cluster}.log",
    conda:
        "../../envs/alignment.yaml"
    shell:
        "julia {params} {input} > {output} 2> {log}"


rule trimal:
    input:
        outdir / "aligned/taper/{probe}/{probe}_{cluster}.fasta",
    output:
        outdir / "aligned/trimal/{probe}/{probe}_{cluster}.fasta",
    params:
        "-gt " + str(trimal_gt),  # discard columns that have more than than 80% of gaps
    log:
        outdir / "logs/aligned/taper/{probe}/{probe}_{cluster}.log",
    conda:
        "../../envs/alignment.yaml"
    shell:
        "trimal {params} -in {input} -out {output} 2> {log}"


def aggregate_clean_alns(wildcards):
    checkpoint_output = checkpoints.separate_multiseq.get(probe=wildcards.probe).output[
        0
    ]
    glob_match = glob_wildcards(
        Path(checkpoint_output) / f"{wildcards.probe}_{{cluster}}.fasta"
    )
    return expand(
        outdir / "aligned/trimal/{probe}/{probe}_{cluster}.fasta",
        probe=wildcards.probe,
        cluster=glob_match.cluster,
    )


rule collect_clean_alignments:
    input:
        aggregate_clean_alns,
    output:
        chkpt=outdir / "logs/aligned/collect_clean_alns/{probe}.chkpt",
    shell:
        "echo {input} | tr '[:space:]' '\n' >> {output.chkpt}"


checkpoint done_cleaning:
    input:
        chkpt=expand(
            outdir / "logs/aligned/collect_clean_alns/{probe}.chkpt", probe=probes_list
        ),
    output:
        done=touch(outdir / "logs/dones/clean_alns.done"),
