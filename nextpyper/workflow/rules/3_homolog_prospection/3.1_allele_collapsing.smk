rule allele_clustering:
    input:
        outdir / "saute/final/merged/{sample}.fasta",
    output:
        outdir / "homolog_prospection/allele_collapsing/{sample}.fasta",
    log:
        outdir / "logs/homolog_prospection/allele_collapsing/{sample}.log",
    params:
        "--id 0.99 --minseqlength 5 --qmask none",
    threads: 4
    conda:
        "../../envs/clustering.yaml"
    shell:
        """
        vsearch --threads {threads} {params} \
            --cluster_fast {input} \
            --centroids {output} \
            2> {log}
        """
