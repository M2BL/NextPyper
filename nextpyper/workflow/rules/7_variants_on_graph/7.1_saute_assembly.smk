rule merge_consensus_probes:
    input:
        expand(outdir / "clustering/consensus/{probes}.fasta", probes=probes_list),
    output:
        outdir / "saute/consensus.fasta",
    shell:
        "cat {input} > {output}"


rule saute_assembly:
    input:
        reads1=outdir / "preprocessed/cleaned/{sample}_R1.fastq.gz",
        reads2=outdir / "preprocessed/cleaned/{sample}_R2.fastq.gz",
        consensus=outdir / "saute/consensus.fasta",
    output:
        all_vars=outdir / "saute/target_assembly/{sample}/all_vars.fasta",
        target_vars=outdir / "saute/target_assembly/{sample}/target_vars.fasta",
        graph=outdir / "saute/target_assembly/{sample}/graph.gfa",
    params:
        "--max_variants 10000 --extend_ends --remove_homopolymer_indels",
    log:
        outdir / "logs/saute/{sample}.log",
    threads: 8
    conda:
        "../../envs/saute.yaml"
    shell:
        "(saute --cores {threads} {params} "
        "--reads {input.reads1},{input.reads2} "
        "--targets {input.consensus} "
        "--gfa {output.graph} "
        "--all_variants {output.all_vars} "
        "--selected_variants {output.target_vars}) > {log} 2>&1 "
