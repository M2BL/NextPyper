def aggregate_sample_per_probe(wildcards):
    probe_inputs = defaultdict(list)
    for sample in sample_list:
        checkpoint_output = checkpoints.homologs_filtering.get(sample=sample).output[0]
        global_match = glob_wildcards(Path(checkpoint_output) / "{probe}.fasta")

        for probe in global_match.probe:
            probe_inputs[probe].append(
                outdir / f"assembled/filtering/filtered_scfs/{sample}/{probe}.fasta"
            )

    return probe_inputs[wildcards.probe]


rule merge_asms:
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
        centroids=outdir / "clustering/clusters/{probe}.fasta",
    log:
        outdir / "logs/clustering/vsearch/{probe}.log",
    params:
        extra="--id 0.95 --iddef 3 --minseqlength 5 --qmask none --strand both",
    threads: 1
    wrapper:
        "v4.3.0/bio/vsearch"


rule vsearch_renaming:
    input:
        outdir / "clustering/clusters/{probe}.fasta",
    output:
        outdir / "clustering/consensus/{probe}.fasta",
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """awk '/^>/ {{header=$0; sub(/^>/, "", header); 
        split(header, parts, "_"); 
        printf ">"; 
        for(i=1; i<=length(parts)-4; i++) 
        {{printf "%s%s", (i>1?"_":""), parts[i]}}; 
        print ""}} !/^>/ {{print}}' {input} > {output}
        """
