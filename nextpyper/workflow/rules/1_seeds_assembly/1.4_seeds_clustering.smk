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
        cluster_fast=outdir / "clustering/sample_merged_input/{probe}.fasta",
    output:
        uc=outdir / "clustering/cluster_tables/{probe}.tsv",
    log:
        outdir / "logs/clustering/vsearch/{probe}.log",
    params:
        extra="--id 0.95 --iddef 3 --minseqlength 5 --qmask none --strand both",
    threads: 4
    wrapper:
        "v4.3.0/bio/vsearch"


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
        spades_folders=expand(outdir / "assembled/spades/{sample}", sample=sample_list),
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
        min_sister_freq=min_sister_sample_freq,
        pattern=pattern,
        is_multi=multi_probes,
        interseeds_use=interseeds_use,
        cov_by_mapping=lookup("seeds/cov_by_mapping", within=pipeline),
        heuristic_params=saute_heuristic_params,
    script:
        "../../../src/seeds_collection.py"
