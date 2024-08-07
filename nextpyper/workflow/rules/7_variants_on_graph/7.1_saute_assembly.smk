targets.append(
    expand(outdir / "saute/target_assembly/{samples}/graph.gfa", samples=sample_list)
)


def gather_consensus(wildcards):
    checkpoint_output = checkpoints.done_hmms.get(**wildcards).output[0]
    glob_match = glob_wildcards(
        outdir / "HMM/consensus/{probe}/{probe2}_{cluster}.fasta"
    )

    return expand(
        outdir / "HMM/consensus/{probe}/{probe}_{cluster}.fasta",
        zip,
        probe=glob_match.probe,
        cluster=glob_match.cluster,
    )


rule merge_consensus_probes:
    input:
        direc=outdir / "logs/dones/hmms.done",
        consensus=gather_consensus,
    output:
        outdir / "saute/consensus.fasta",
    shell:
        "cat {input.consensus} > {output}"


rule saute_assembly:
    input:
        reads1=outdir / "trimmed/{sample}_R1.fastq",
        reads2=outdir / "trimmed/{sample}_R2.fastq",
        consensus=outdir / "saute/consensus.fasta",
    output:
        all_vars=outdir / "saute/target_assembly/{sample}/all_vars.fasta",
        target_vars=outdir / "saute/target_assembly/{sample}/target_vars.fasta",
        graph=outdir / "saute/target_assembly/{sample}/graph.gfa",
    params:
        "--max_variants 10000 ",
    log:
        outdir / "logs/saute/{sample}.log",
    threads: 8
    conda:
        "../../envs/saute.yaml"
    shell:
        "(saute --cores {threads} {params} "
        "--reads {input.reads1},{input.reads2} "
        "--targets {input.consensus} "
        "--gfa {output.graph} "
        "--all_variants {output.all_vars} "
        "--selected_variants {output.target_vars}) > {log} 2>&1 "
