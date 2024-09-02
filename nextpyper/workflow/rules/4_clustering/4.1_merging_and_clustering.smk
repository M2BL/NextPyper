graph_dir = outdir / Path("assembled/prefixed/component_seqs/")


targets.append(expand(outdir / "clustering/clusters/{probes}", probes=probes_list))
# ToDo: Test:
# What happens if not all the probes are present in at least one sample?
# This would certainly happen in real life cases.


def union_probes(wildcards):
    glob_match = glob_wildcards(graph_dir / f"{{sample}}/{wildcards.probe}.fasta")

    return expand(
        graph_dir / f"{{sample}}/{wildcards.probe}.fasta",
        sample=glob_match.sample,
    )


rule merge_asms:
    input:
        direc=outdir / "logs/dones/prefixing.done",
        probes=union_probes,
    output:
        outfile=outdir / "clustering/sample_merged_input/{probe}.fasta",
    shell:
        "cat {input.probes} > {output}"


checkpoint clustering:
    input:
        probes=outdir / "translated_probes/split_probes/{probe}.fasta",
        contigs=outdir / "clustering/sample_merged_input/{probe}.fasta",
    output:
        directory(outdir / "clustering/clusters/{probe}"),
    log:
        outdir / "logs/clustering/{probe}.log",
    conda:
        "../../envs/clustering.yaml"
    script:
        "../../../src/contig_cluster.py"


rule vsearch_clustering:
    input:
        cluster_fast=outdir / "clustering/sample_merged_input/{probe}.fasta",
    output:
        msaout=outdir / "clustering/clusters/{probe}.fasta",
    log:
        outdir / "logs/clustering/{probe}.log",
    params:
        extra="--id 0.95 --minseqlength 5",
    threads: 1
    wrapper:
        "v4.3.0/bio/vsearch"
