from vsearch import get_vsearch_kmer_consensus


def union_probes(wildcards):
    if not multi_probes:
        glob_match = glob_wildcards(
            outdir / f"assembled/split_components/{{sample}}/{wildcards.probe}.fasta"
        )
        return expand(
            outdir / f"assembled/split_components/{{sample}}/{wildcards.probe}.fasta",
            sample=glob_match.sample,
        )
    else:
        glob_match = glob_wildcards(
            outdir
            / f"assembled/split_components/{{sample}}/{wildcards.probe}_{{cluster}}.fasta"
        )

        return expand(
            outdir
            / f"assembled/split_components/{{sample}}/{wildcards.probe}_{{cluster}}.fasta",
            zip,
            sample=glob_match.sample,
            cluster=glob_match.cluster,
        )


rule merge_asms:
    input:
        direc=outdir / "logs/dones/splitting.done",
        probes=union_probes,
    output:
        outfile=outdir / "clustering/sample_merged_input/{probe}.fasta",
    shell:
        """
        for file in {input.probes}; do 
            cat $file 
        done > {output}
        touch {output}
        """


rule vsearch_clustering:
    input:
        cluster_fast=outdir / "clustering/sample_merged_input/{probe}.fasta",
    output:
        msaout=outdir / "clustering/clusters/{probe}.fasta",
    log:
        outdir / "logs/clustering/vsearch/{probe}.log",
    params:
        extra="--id 0.95 --minseqlength 5 --qmask none",
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
