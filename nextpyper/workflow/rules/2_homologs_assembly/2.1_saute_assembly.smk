def low_cov_params(wildcards, input):
    if (med_cov := float(Path(input.cov).read_text())) > 20:
        return ""
    elif med_cov < 10:
        k = 49
    else:
        with open(input.json) as file:
            summary = json.load(file)
            k_mid = int(summary["summary"]["after_filtering"]["read2_mean_length"]) // 2
            k = k_mid + 1 if k_mid % 2 == 0 else k_mid

    return f"--kmer {k} --secondary_kmer 21 --secondary_kmer_threshold 5"


rule saute_assembly:
    input:
        reads1=outdir / "preprocessed/cleaned/{sample}_R1.fastq.gz",
        reads2=outdir / "preprocessed/cleaned/{sample}_R2.fastq.gz",
        seeds=outdir / "saute/seeds/{sample}.fasta",
        cov=outdir / "logs/clustering/seed_collection/{sample}.cov",
        json=outdir / "logs/preprocessing/fastp/{sample}.json",
    output:
        all_vars=outdir / "saute/target_assembly/{sample}/all_vars.fasta",
        target_vars=outdir / "saute/target_assembly/{sample}/target_vars.fasta",
        graph=outdir / "saute/target_assembly/{sample}/graph.gfa",
    params:
        glob="--max_variants 10 --extend_ends --remove_homopolymer_indels ",
        target_cov=saute_target_cov,
        cov=low_cov_params,
    log:
        outdir / "logs/saute/{sample}.log",
    threads: 8
    conda:
        "../../envs/saute.yaml"
    shell:
        "(saute --cores {threads} {params.glob} {params.cov} "
        "--target_coverage {params.target_cov} "
        "--reads {input.reads1},{input.reads2} "
        "--targets {input.seeds} "
        "--gfa {output.graph} "
        "--all_variants {output.all_vars} "
        "--selected_variants {output.target_vars}) > {log} 2>&1 "


rule fix_homologs_header:
    input:
        outdir / "saute/target_assembly/{sample}/all_vars.fasta",
    output:
        outdir / "saute/target_assembly/{sample}/fixed_vars.fasta",
    params:
        pattern=SAUTE_PRE_FIX_PAT,
        sample=lambda wildcards: wildcards.sample,
    script:
        "../../../src/fix_headers.py"
