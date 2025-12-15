rule distribute_seeds:
    input:
        expand(
            outdir / "assembled/filtering/filtered_scfs/{sample}.fasta",
            sample=sample_list,
        ),
    output:
        expand(
            outdir / "clustering/sample_merged_input/{probe}.fasta", probe=probes_list
        ),
    log:
        outdir / "logs/clustering/seed_distribution.log",
    params:
        pattern=lambda wildcards: r"-(?P<probe>.*?)_",
        probes=probes_list,
        mode="scfs",
    script:
        "../../../src/multi_seq_probes.py"


rule vsearch_clustering:
    input:
        outdir / "clustering/sample_merged_input/{probe}.fasta",
    output:
        outdir / "clustering/cluster_tables/{probe}.tsv",
    log:
        outdir / "logs/clustering/vsearch/{probe}.log",
    params:
        "--id 0.95 --iddef 3 --minseqlength 5 --qmask none --strand both    ",
    threads: 4
    conda:
        "../../envs/clustering.yaml"
    shell:
        """
        vsearch --threads {threads} {params} \
            --cluster_fast {input} \
            --uc {output} \
            2> {log}
        """


rule seeds_collection:
    input:
        probes=probes_path.resolve(),
        cluster_tables=expand(
            outdir / "clustering/cluster_tables/{probe}.tsv", probe=probes_list
        ),
        samples=expand(
            outdir / "assembled/filtering/filtered_scfs/{sample}.fasta",
            sample=sample_list,
        ),
        spades_graphs=expand(
            outdir / "assembled/spades/{sample}/assembly_graph_with_scaffolds.gfa",
            sample=sample_list,
        ),
        covs=expand(
            outdir / "assembled/filtering/coverage/{sample}.metabat",
            sample=sample_list,
        ),
        read_stats=expand(
            outdir / "logs/preprocessing/fastp/{sample}.json", sample=sample_list
        ),
    output:
        seeds=expand(outdir / "saute/seeds/{sample}.fasta", sample=sample_list),
        saute_params=expand(
            outdir / "logs/saute/kmer_params/{sample}.json", sample=sample_list
        ),
    params:
        min_sister_freq=lookup("seeds/min_sister_sample_freq", within=pipeline),
        pattern=probe_pattern,
        is_multi=multi_probes,
        interseeds_use=interseeds_use,
        cov_by_mapping=lookup("seeds/cov_by_mapping", within=pipeline),
        heuristic_params=lookup("saute/heuristic", within=pipeline),
    script:
        "../../../src/seeds_collection.py"
