pathvars:
    stage="1_assembly/13_clustering",
    results="<workdir>/<stage>",
    inter="1_assembly/12_filtering/122_filt_seeds/1221_filtered_scfs"

    
rule distribute_seeds:
    input:
        expand("<workdir>/<inter>/{sample}.fasta", sample=sample_list),
    output:
        expand("<results>/130_sample_merged_input/{probe}.fasta", probe=probes_list),
    log:
        "<logs>/<stage>/seed_distribution.log",
    params:
        pattern=lambda wildcards: r"-(?P<probe>.*?)_EDGE_",
        probes=probes_list,
        mode="scfs",
    script:
        "../../../src/multi_seq_probes.py"


rule vsearch_clustering:
    input:
        "<results>/130_sample_merged_input/{probe}.fasta"
    output:
        "<results>/131_cluster_tables/{probe}.tsv",
    log:
        "<logs>/<stage>/131_vsearch/{probe}.log",
    conda:
        "../../envs/clustering.yaml"
    threads: 4
    params:
        "--id 0.95 --iddef 3 --minseqlength 5 --qmask none --strand both",
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
        cluster_tables=expand("<results>/131_cluster_tables/{probe}.tsv", probe=probes_list),
        samples=expand("<workdir>/<inter>/{sample}.fasta", sample=sample_list),
        spades_graphs=expand(
            "<workdir>/1_assembly/10_raw_assembly/100_spades/{sample}/assembly_graph_with_scaffolds.gfa",
            sample=sample_list,
        ),
        covs=expand(
            "<workdir>/1_assembly/12_filtering/121_coverage/{sample}.metabat",
            sample=sample_list,
        ),
        read_stats=expand(
            "<logs>/0_preprocessing/01_trimming_fastp/{sample}.json", sample=sample_list
        ),
    output:
        seeds=expand("<workdir>/2_saute/21_seeds/{sample}.fasta", sample=sample_list),
        saute_params=expand("<workdir>/2_saute/20_params_and_stats/200_kmer_params/{sample}.json", sample=sample_list),
    log:
        "<logs>/2_saute/seeds_collection.log",
    params:
        pattern=probe_pattern,
        is_multi=multi_probes,
        interseeds_use=interseeds_use,
        min_sister_freq=lookup("seeds/min_sister_sample_freq", within=pipeline),
        cov_by_mapping=lookup("seeds/cov_by_mapping", within=pipeline),
        heuristic_params=lookup("saute/heuristic", within=pipeline),
    script:
        "../../../src/seeds_collection.py"
