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
        extra="--extend_ends --remove_homopolymer_indels ",
        kmer_threshold=lookup(
            "saute/assembly/secondary_kmer_threshold", within=pipeline
        ),
        max_var=lookup("saute/assembly/max_variants", within=pipeline),
        target_cov=lookup("saute/assembly/target_cov", within=pipeline),
        kmers=saute_kmer,
    log:
        outdir / "logs/saute/assembly/{sample}.log",
    threads: 8
    conda:
        "../../envs/saute.yaml"
    shell:
        """saute --cores {threads} {params.extra} {params.kmers} \
        --target_coverage {params.target_cov} \
        --max_variants {params.max_var} \
        --secondary_kmer_threshold {params.kmer_threshold} \
        --reads {input.reads1},{input.reads2} \
        --targets {input.seeds} \
        --gfa {output.graph} \
        --all_variants {output.all_vars} \
        --selected_variants {output.target_vars} > {log} 2>&1
        """


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
