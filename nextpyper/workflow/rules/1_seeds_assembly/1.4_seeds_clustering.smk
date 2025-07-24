def aggregate_sample_per_probe(wildcards):
    probe_inputs = defaultdict(list)
    for sample in sample_list:
        checkpoint_output = checkpoints.seeds_filtering.get(sample=sample).output[0]
        global_match = glob_wildcards(Path(checkpoint_output) / "{probe}.fasta")

        for probe in global_match.probe:
            probe_inputs[probe].append(
                outdir / f"assembled/filtering/filtered_scfs/{sample}/{probe}.fasta"
            )

    return probe_inputs[wildcards.probe]


rule merge_samples_seeds:
    input:
        probe=aggregate_sample_per_probe,
        chkpt=outdir / "logs/dones/splitting.done",
    output:
        outfile=outdir / "clustering/sample_merged_input/{probe}.fasta",
    shell:
        """
        for file in {input.probe}; do 
            cat $file 
        done > {output}
        touch {output}
        """


rule vsearch_clustering:
    input:
        cluster_fast=outdir / "clustering/sample_merged_input/{probe}.fasta",
    output:
        centroids=outdir / "clustering/centroids/{probe}.fasta",
        msaout=outdir / "clustering/msa/{probe}.fasta",
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
            outdir / "assembled/filtering/filtered_scfs/{sample}",
            sample=sample_list,
        ),
        spades_folders=expand(outdir / "assembled/spades/{sample}", sample=sample_list),
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
        heuristic_params=saute_heuristic_params,
    script:
        "../../../src/seeds_collection.py"
