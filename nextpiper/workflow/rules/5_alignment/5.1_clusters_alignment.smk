from itertools import chain

targets.append(outdir / "logs/dones/cluster_alns.done")


def nrecs(fasta):
    return sum(1 for _ in SeqIO.parse(fasta, "fasta"))


checkpoint separate_multiseq:
    input:
        outdir / "clustering/clusters/{probe}",
    output:
        directory(outdir / "aligned/aln_inputs/{probe}"),
    run:
        for file in Path(input[0]).glob("*.fasta"):
            if nrecs(file) > 1:
                Path(output[0]).mkdir(parents=True, exist_ok=True)
                outfile = Path(output[0]) / file.name
                outfile.symlink_to(file.resolve())


def get_aln_input(wildcards):
    checkpoint_output = checkpoints.separate_multiseq.get(**wildcards).output[0]
    return Path(checkpoint_output) / f"{wildcards.probe}_{wildcards.cluster}.fasta"


rule mafft:
    input:
        clusters=get_aln_input,
    output:
        alns=outdir / "aligned/cluster_alns/{probe}/{probe}_{cluster}.fasta",
    params:
        "--auto --adjustdirection",
    log:
        outdir / "logs/aligned/cluster_alns/{probe}/{probe}_{cluster}.log",
    conda:
        "../../envs/alignment.yaml"
    shell:
        "mafft {params} {input.clusters} > {output.alns} 2> {log}"


def aggregate_alns(wildcards):
    checkpoint_output = checkpoints.separate_multiseq.get(probe=wildcards.probe).output[
        0
    ]
    glob_match = glob_wildcards(
        Path(checkpoint_output) / f"{wildcards.probe}_{{cluster}}.fasta"
    )
    return expand(
        outdir / "aligned/cluster_alns/{probe}/{probe}_{cluster}.fasta",
        probe=wildcards.probe,
        cluster=glob_match.cluster,
    )


rule collect_alignments:
    input:
        aggregate_alns,
    output:
        chkpt=outdir / "logs/aligned/collect_alns/{probe}.chkpt",
    shell:
        "echo {input} | tr '[:space:]' '\n' >> {output.chkpt}"


checkpoint done_alns:
    input:
        chkpt=expand(
            outdir / "logs/aligned/collect_alns/{probe}.chkpt", probe=probes_list
        ),
    output:
        done=touch(outdir / "logs/dones/cluster_alns.done"),
