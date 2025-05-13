use rule vsearch_clustering as allele_clustering with:
    input:
        cluster_fast=outdir
        / "homolog_prospection/homologs_filtering/filtered_scfs/{sample}.fasta",
    output:
        centroids=outdir / "homolog_prospection/allele_collapsing/{sample}.fasta",
    log:
        outdir / "logs/homolog_prospection/allele_collapsing/{sample}.log",
    params:
        extra="--id 0.99 --minseqlength 5 --qmask none",
