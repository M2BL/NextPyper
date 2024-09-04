from vsearch import get_vsearch_kmer_consensus

targets.append(
    expand(outdir / "clustering/consensus/{probes}.fasta", probes=probes_list)
)
# ToDo: Test:
# What happens if not all the probes are present in at least one sample?
# This would certainly happen in real life cases.


def union_probes(wildcards):
    glob_match = glob_wildcards(
        outdir / f"assembled/split_components/{{sample}}/{wildcards.probe}.fasta"
    )

    return expand(
        outdir / f"assembled/split_components/{{sample}}/{wildcards.probe}.fasta",
        sample=glob_match.sample,
    )


rule merge_asms:
    input:
        direc=outdir / "logs/dones/splitting.done",
        probes=union_probes,
    output:
        outfile=outdir / "clustering/sample_merged_input/{probe}.fasta",
    shell:
        "cat {input.probes} > {output}"


rule vsearch_clustering:
    input:
        cluster_fast=outdir / "clustering/sample_merged_input/{probe}.fasta",
    output:
        msaout=outdir / "clustering/clusters/{probe}.fasta",
    log:
        outdir / "logs/clustering/vsearch/{probe}.log",
    params:
        extra="--id 0.95 --minseqlength 5",
    threads: 1
    wrapper:
        "v4.3.0/bio/vsearch"


rule vsearch_consensus_parsing:
    input:
        outdir / "clustering/clusters/{probe}.fasta",
    output:
        outdir / "clustering/consensus/{probe}.fasta",
    log:
        outdir / "logs/clustering/consensus/{probe}.log",
    run:
        with open(log[0], "w") as outlog:
            sys.stdout = sys.stderr = outlog
            recs = get_vsearch_kmer_consensus(Path(input[0]), "SPAdes")
            SeqIO.write(recs, Path(output[0]), "fasta")
