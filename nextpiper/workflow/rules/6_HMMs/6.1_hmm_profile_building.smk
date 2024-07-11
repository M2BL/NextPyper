from hmm_build import hmm_build, hmm_consensus

targets.append(outdir / "logs/dones/hmms.done")


def is_multiseq(wildcards) -> bool:
    return 1 < nrecs(
        outdir
        / f"clustering/clusters/{wildcards.probe}/{wildcards.probe}_{wildcards.cluster}.fasta"
    )


rule build_hmms:
    input:
        branch(
            is_multiseq,
            then=outdir / "aligned/trimal/{probe}/{probe}_{cluster}.fasta",
            otherwise=outdir / "clustering/clusters/{probe}/{probe}_{cluster}.fasta",
        ),
    output:
        prof=outdir / "HMM/profiles/{probe}/{probe}_{cluster}.hmm",
        consensus=outdir / "HMM/consensus/{probe}/{probe}_{cluster}.fasta",
    log:
        outdir / "logs/HMM/build/{probe}/{probe}_{cluster}.log",
    run:
        hmm_build(Path(input[0]), Path(output.prof), "amino")
        hmm_consensus(Path(output.prof), Path(output.consensus))


def aggregate_hmms(wildcards):
    checkpoint_output = checkpoints.clustering.get(probe=wildcards.probe).output[0]
    glob_match = glob_wildcards(
        Path(checkpoint_output) / f"{wildcards.probe}_{{cluster}}.fasta"
    )

    return expand(
        outdir / "HMM/profiles/{probe}/{probe}_{cluster}.hmm",
        probe=wildcards.probe,
        cluster=glob_match.cluster,
    )


rule collect_hmms:
    input:
        aggregate_hmms,
    output:
        chkpt=outdir / "logs/HMM/collect_hmms/{probe}.chkpt",
    shell:
        "echo {input} | tr '[:space:]' '\n' >> {output.chkpt}"


checkpoint done_hmms:
    input:
        chkpt=expand(outdir / "logs/HMM/collect_hmms/{probe}.chkpt", probe=probes_list),
    output:
        done=touch(outdir / "logs/dones/hmms.done"),
