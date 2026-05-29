pathvars:
    stage="3_homolog_prospection/30_allele_collapsing",
    results="<workdir>/<stage>",

rule allele_clustering:
    input:
        "<workdir>/2_saute/24_postprocessing/242_reheading/{sample}.fasta",
    output:
        "<results>/{sample}.fasta",
    log:
        "<logs>/<stage>/{sample}.log",
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
