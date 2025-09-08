def read_kmer_params(kmer_params: Path) -> tuple[int, int]:
    kmer_params = json.loads(kmer_params.read_text())
    k1 = int(kmer_params["k1"])
    k2 = int(kmer_params["k2"])
    return k1, k2


def saute_kmer(wildcards, input):
    """Set the Kmer sizes for saute assembly"""

    k1, k2 = read_kmer_params(Path(input.kmer_params))
    return f"--kmer {k1} --secondary_kmer {k2}"


def saute_kmer_expl(wildcards, input):
    """Set the Kmer sizes for saute reassembly rescaling the primary kmer"""

    k1rescale = lookup("saute/reassembly/k1_rescaling", within=pipeline)
    k1, k2 = read_kmer_params(Path(input.kmer_params))

    ## For reassembly, scale the primary kmer
    if k1rescale:
        k1 = int(k1 * k1rescale)
        k1 = k1 - 1 if k1 % 2 == 0 else k1

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
        max_var=lookup("saute/max_variants", within=pipeline),
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


checkpoint split_saute_assembly:
    input:
        outdir / "saute/target_assembly/{sample}/target_vars.fasta",
    output:
        normal=outdir / "saute/target_assembly/{sample}/normal_vars.fasta",
        expl=outdir / "saute/target_assembly/{sample}/expl_vars.fasta",
    params:
        pattern=TARGET_COLLAPSE_PAT,
        explosive_limit=lookup("saute/reassembly/explosive_limit", within=pipeline),
    log:
        outdir / "logs/saute/reassembly/split/{sample}.log",
    script:
        "../../../src/var_asm_parser.py"


rule collect_explosive_reads:
    input:
        scfs=outdir / "saute/target_assembly/{sample}/expl_vars.fasta",
        reads1=outdir / "preprocessed/cleaned/{sample}_R1.fastq.gz",
        reads2=outdir / "preprocessed/cleaned/{sample}_R2.fastq.gz",
    output:
        reads1=outdir / "saute/expl_assembly/{sample}/expl_R1.fastq.gz",
        reads2=outdir / "saute/expl_assembly/{sample}/expl_R2.fastq.gz",
    log:
        outdir / "logs/saute/reassembly/read_collection/{sample}.log",
    threads: 4
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """minimap2 -t {threads} -ax sr {input.scfs} {input.reads1} {input.reads2} 2> {log} | \
        samtools view -@ {threads} -uhf 2 2>> {log} | samtools sort -n -u -@ {threads} 2>> {log} | \
        samtools fastq -1 {output.reads1} -2 {output.reads2} - 2>> {log}
        """


rule collect_explosive_seeds:
    input:
        seeds=outdir / "saute/seeds/{sample}.fasta",
        expl=outdir / "saute/target_assembly/{sample}/expl_vars.fasta",
    output:
        seeds=outdir / "saute/expl_assembly/{sample}/expl_seeds.fasta",
        seqids=temp(outdir / "saute/expl_assembly/{sample}/seqids.txt"),
    log:
        outdir / "logs/saute/reassembly/seed_collection/{sample}.log",
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """awk '/>/{{match($0,/-(.*?)_EDGE/, m); print "-"m[1]"_"}}' {input.expl} | sort | uniq > {output.seqids}
        seqkit grep -rf {output.seqids} {input.seeds} > {output.seeds}
        """


rule explosive_reassembly:
    input:
        reads1=outdir / "saute/expl_assembly/{sample}/expl_R1.fastq.gz",
        reads2=outdir / "saute/expl_assembly/{sample}/expl_R2.fastq.gz",
        seeds=outdir / "saute/expl_assembly/{sample}/expl_seeds.fasta",
        kmer_params=outdir / "logs/saute/kmer_params/{sample}.json",
    output:
        all_vars=outdir / "saute/expl_assembly/{sample}/all_vars.fasta",
        target_vars=outdir / "saute/expl_assembly/{sample}/target_vars.fasta",
        graph=outdir / "saute/expl_assembly/{sample}/graph.gfa",
    params:
        extra="--extend_ends --remove_homopolymer_indels ",
        kmer_threshold=lookup(
            "saute/reassembly/secondary_kmer_threshold", within=pipeline
        ),
        max_var=lookup("saute/max_variants", within=pipeline),
        target_cov=lookup("saute/reassembly/target_cov", within=pipeline),
        k1rescale=lookup("saute/reassembly/k1_rescaling", within=pipeline),
        kmers=saute_kmer_expl,
    log:
        outdir / "logs/saute/reassembly/assembly/{sample}.log",
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
        --selected_variants {output.target_vars} > {log} 2>&1 || \
        touch {output.all_vars} {output.graph} {output.target_vars}
        """


checkpoint collapse_alleles_explosive:
    input:
        outdir / "saute/expl_assembly/{sample}/target_vars.fasta",
    output:
        normal=outdir / "saute/expl_assembly/{sample}/collapsed_vars.fasta",
    params:
        pattern=TARGET_COLLAPSE_PAT,
        empty_ok=True,
    log:
        outdir / "logs/saute/reassembly/expl_collapse_alleles/{sample}.log",
    script:
        "../../../src/var_asm_parser.py"


rule normal_vars_check:
    input:
        outdir / "saute/target_assembly/{sample}/expl_vars.fasta",
    output:
        touch(outdir / "saute/expl_assembly/{sample}/all_normal.chkp"),


# All probes are normal, no need to do reassembly.
def all_normal(wildcards):
    out_expl = checkpoints.split_saute_assembly.get(sample=wildcards.sample).output.expl
    return Path(out_expl).stat().st_size == 0


# Reassembly yield nothing, so take back the original results.
def empty_explosive_asm(wildcards):
    out_expl = checkpoints.collapse_alleles_explosive.get(
        sample=wildcards.sample
    ).output.normal
    return Path(out_expl).stat().st_size == 0


rule collect_saute_assemblies:
    input:
        normal=outdir / "saute/target_assembly/{sample}/normal_vars.fasta",
        expl=branch(
            all_normal,
            then=outdir / "saute/expl_assembly/{sample}/all_normal.chkp",
            otherwise=branch(
                not reasm or empty_explosive_asm,
                then=outdir / "saute/target_assembly/{sample}/expl_vars.fasta",
                otherwise=outdir / "saute/expl_assembly/{sample}/collapsed_vars.fasta",
            ),
        ),
    output:
        outdir / "saute/final/collected/{sample}.fasta",
    shell:
        "cat {input.normal} {input.expl} > {output}"


rule collapse_variants:
    input:
        outdir / "saute/final/collected/{sample}.fasta",
    output:
        normal=outdir / "saute/final/collapsed/{sample}.fasta",
    params:
        pattern=TARGET_COLLAPSE_PAT,
        collapse_vars=lookup("saute/collapse_vars", within=pipeline),
        max_var=lookup("saute/max_variants", within=pipeline),
    log:
        outdir / "logs/saute/variant_collapsing/{sample}.log",
    script:
        "../../../src/var_asm_parser.py"


rule fix_homologs_header:
    input:
        branch(
            lookup("saute/collapse_vars", within=pipeline),
            then=outdir / "saute/final/collapsed/{sample}.fasta",
            otherwise=outdir / "saute/final/collected/{sample}.fasta",
        ),
    output:
        outdir / "saute/final/merged/{sample}.fasta",
    params:
        pattern=SAUTE_PRE_FIX_PAT,
        sample=lambda wildcards: wildcards.sample,
    script:
        "../../../src/fix_headers.py"
