from vsearch import get_vsearch_kmer_consensus


rule allele_clustering:
    input:
        cluster_fast=outdir
        / "homolog_prospection/blast_filtering/filtered_scfs/{sample}.fasta",
    output:
        msaout=outdir / "homolog_prospection/allele_collapsing/vsearch/{sample}.fasta",
    log:
        outdir / "logs/homolog_prospection/allele_collapsing/vsearch/{sample}.log",
    params:
        extra="--id 0.99 --minseqlength 5 --qmask none",
    threads: 1
    wrapper:
        "v4.3.0/bio/vsearch"


rule allele_parsing:
    input:
        outdir / "homolog_prospection/allele_collapsing/vsearch/{sample}.fasta",
    output:
        outdir / "homolog_prospection/allele_collapsing/consensus/{sample}.fasta",
    log:
        outdir / "logs/homolog_prospection/allele_collapsing/consensus/{sample}.log",
    run:
        with open(log[0], "w") as outlog:
            sys.stdout = sys.stderr = outlog
            recs = get_vsearch_kmer_consensus(Path(input[0]), "SAUTE")
            SeqIO.write(recs, Path(output[0]), "fasta")
