def saute_kmer_params(wildcards, input):
    # Read the median kmer cov in spades scaffolds
    med_cov = float(Path(input.cov).read_text())
    # Read the read lenght of reads
    with open(input.json) as file:
        summary = json.load(file)
        L = int(summary["summary"]["after_filtering"]["read2_mean_length"])

    # Read the Kmer size used by spades
    spades_folder = outdir / f"assembled/spades/{wildcards.sample}"
    kspades = int(one(spades_folder.glob("K*/final_contigs.paths")).parent.name[1:])

    # Compute median depth observed
    read_med_cov = med_cov * L / (L - kspades + 1)

    # Compute kmer sizes for the given targets
    k2nd = int(L * (1 - (secondary_target_depth / read_med_cov)) + 1)
    k1st = int(L * (1 - (primary_target_depth / read_med_cov)) + 1)

    # Adjust kmer to ranges kmer: 1st [0.5 - 0.8]L, 2nd [0.15 - 0.3]L
    k2nd_min, k2nd_max = int(L * 0.15), int(L * 0.3)
    k1st_min, k1st_max = int(L * 0.5), int(L * 0.8)

    if k2nd < k2nd_min:
        k2nd = k2nd_min
    elif k2nd > k2nd_max:
        k2nd = k2nd_max

    if k1st < k1st_min:
        k1st = k1st_min
    elif k1st > k1st_max:
        k1st = k1st_max

    # Cap secondary kmer on 21 and ensure odd kmers
    k2nd = 21 if k2nd < 21 else k2nd
    k2nd = k2nd + 1 if k2nd % 2 == 0 else k2nd
    k1st = k1st + 1 if k1st % 2 == 0 else k1st

    return f"--kmer {k1st} --secondary_kmer {k2nd}"


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
        glob="--max_variants 10 --extend_ends --remove_homopolymer_indels --secondary_kmer_threshold 3 ",
        target_cov=saute_target_cov,
        cov=saute_kmer_params,
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
