def saute_kmer(wildcards, input):
    with open(input.kmer_params) as file:
        kmer_params = json.load(file)
        k1 = int(kmer_params["k1"])
        k2 = int(kmer_params["k2"])

    return f"--kmer {k1} --secondary_kmer {k2}"


rule saute_assembly:
    input:
        reads1=outdir / "preprocessed/cleaned/{sample}_R1.fastq.gz",
        reads2=outdir / "preprocessed/cleaned/{sample}_R2.fastq.gz",
        seeds=outdir / "saute/seeds/{sample}.fasta",
        kmer_params=outdir / "logs/saute/kmer_params/{sample}.json",
    output:
        all_vars=outdir / "saute/target_assembly/{sample}/all_vars.fasta",
        target_vars=outdir / "saute/target_assembly/{sample}/target_vars.fasta",
        graph=outdir / "saute/target_assembly/{sample}/graph.gfa",
    params:
        glob="--max_variants 10 --extend_ends --remove_homopolymer_indels --secondary_kmer_threshold 3 ",
        target_cov=saute_target_cov,
        kmers=saute_kmer,
    log:
        outdir / "logs/saute/assembly/{sample}.log",
    threads: 8
    conda:
        "../../envs/saute.yaml"
    shell:
        "(saute --cores {threads} {params.glob} {params.kmers} "
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
