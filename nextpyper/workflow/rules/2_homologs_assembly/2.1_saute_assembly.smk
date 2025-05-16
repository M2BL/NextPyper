rule merge_consensus_probes:
    input:
        expand(outdir / "clustering/centroids/{probes}.fasta", probes=probes_list),
    output:
        outdir / "saute/consensus.fasta",
    shell:
        "cat {input} > {output}"


def aggregate_split(wildcards):
    chkpt_out = checkpoints.seeds_filtering.get(sample=wildcards.sample).output[0]
    return collect(
        outdir
        / f"assembled/filtering/filtered_scfs/{wildcards.sample}/{{probe}}.fasta",
        probe=glob_wildcards(Path(chkpt_out) / "{probe}.fasta").probe,
    )


rule collect_sample_seeds:
    input:
        intra=aggregate_split,
        inter=outdir / "saute/consensus.fasta",
    output:
        outdir / "saute/seeds/{sample}.fasta",
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """
        cat {input.intra} > {output}
        seqkit grep -vrnp {wildcards.sample} {input.inter} >> {output}
        """


rule saute_assembly:
    input:
        reads1=outdir / "preprocessed/trimmed/{sample}_R1.fastq.gz",
        reads2=outdir / "preprocessed/trimmed/{sample}_R2.fastq.gz",
        seeds=outdir / "saute/seeds/{sample}.fasta",
    output:
        all_vars=outdir / "saute/target_assembly/{sample}/all_vars.fasta",
        target_vars=outdir / "saute/target_assembly/{sample}/target_vars.fasta",
        graph=outdir / "saute/target_assembly/{sample}/graph.gfa",
    params:
        "--max_variants 10 --extend_ends --remove_homopolymer_indels",
    log:
        outdir / "logs/saute/{sample}.log",
    threads: 8
    conda:
        "../../envs/saute.yaml"
    shell:
        "(saute --cores {threads} {params} "
        "--reads {input.reads1},{input.reads2} "
        "--targets {input.seeds} "
        "--gfa {output.graph} "
        "--all_variants {output.all_vars} "
        "--selected_variants {output.target_vars}) > {log} 2>&1 "


rule fix_homologs_header:
    input:
        outdir / "saute/target_assembly/{sample}/target_vars.fasta",
    output:
        outdir / "saute/target_assembly/{sample}/fixed_vars.fasta",
    params:
        pattern=SAUTE_PRE_FIX_PAT,
        sample=lambda wildcards: wildcards.sample,
    script:
        "../../../src/fix_headers.py"
