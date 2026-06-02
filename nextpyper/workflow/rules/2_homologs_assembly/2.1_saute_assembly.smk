pathvars:
    stage="2_saute",
    results="<workdir>/<stage>",
    

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
        reads1="<workdir>/0_preprocessing/02_cleaning/{sample}_R1.fastq.gz",
        reads2="<workdir>/0_preprocessing/02_cleaning/{sample}_R2.fastq.gz",
        seeds="<results>/21_seeds/{sample}.fasta",
        kmer_params="<results>/20_params_and_stats/200_kmer_params/{sample}.json",
    output:
        # all_vars=outdir / "saute/target_assembly/{sample}/all_vars.fasta",
        target_vars="<results>/22_homologs_assembly/{sample}/target_vars.fasta",
        graph="<results>/22_homologs_assembly/{sample}/graph.gfa",
    params:
        extra="--extend_ends --remove_homopolymer_indels ",
        kmer_threshold=lookup(
            "saute/assembly/secondary_kmer_threshold", within=pipeline
        ),
        max_var=lookup(query="sample=='{sample}'", cols="homologs", within=sample_table),
        target_cov=lookup("saute/assembly/target_cov", within=pipeline),
        kmers=saute_kmer,
    log:
        "<logs>/<stage>/22_homologs_assembly_saute/{sample}.log",
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
        --selected_variants {output.target_vars} > {log} 2>&1
        """


checkpoint split_saute_assembly:
    input:
        target_vars="<results>/22_homologs_assembly/{sample}/target_vars.fasta",
    output:
        normal="<results>/22_homologs_assembly/{sample}/normal_vars.fasta",
        tribbles="<results>/22_homologs_assembly/{sample}/tribble_vars.fasta",
    params:
        mode="split",
        pattern=TARGET_COLLAPSE_PAT,
        max_vars=lookup(
            query="sample=='{sample}'", cols="homologs", within=sample_table
        ),
    log:
        "<logs>/<stage>/23_tribbles_reassembly/230_split/{sample}.log",
    script:
        "../../../src/var_asm_parser.py"


rule collect_tribble_reads:
    input:
        scfs="<results>/22_homologs_assembly/{sample}/tribble_vars.fasta",
        reads1="<workdir>/0_preprocessing/02_cleaning/{sample}_R1.fastq.gz",
        reads2="<workdir>/0_preprocessing/02_cleaning/{sample}_R2.fastq.gz",
    output:
        reads1="<results>/23_tribbles_reassembly/{sample}/trib_reads_R1.fastq.gz",
        reads2="<results>/23_tribbles_reassembly/{sample}/trib_reads_R2.fastq.gz",
    log:
        "<logs>/<stage>/23_tribbles_reassembly/232_read_collection/{sample}.log",
    threads: 4
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """minimap2 -t {threads} -ax sr {input.scfs} {input.reads1} {input.reads2} 2> {log} | \
        samtools view -@ {threads} -uhf 2 2>> {log} | samtools sort -n -u -@ {threads} 2>> {log} | \
        samtools fastq -1 {output.reads1} -2 {output.reads2} - 2>> {log}
        """


rule collect_tribble_seeds:
    input:
        seeds="<results>/21_seeds/{sample}.fasta",
        tribbles="<results>/22_homologs_assembly/{sample}/tribble_vars.fasta",
    output:
        seeds="<results>/23_tribbles_reassembly/{sample}/tribbles_seeds.fasta",
        seqids=temp("<results>/23_tribbles_reassembly/{sample}/seqids.txt"),
    log:
        "<logs>/<stage>/23_tribbles_reassembly/231_seed_collection/{sample}.log",
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """awk '/>/{{match($0,/-(.*?)_EDGE/, m); print "-"m[1]"_"}}' {input.tribbles} | sort | uniq > {output.seqids}
        seqkit grep -rf {output.seqids} {input.seeds} > {output.seeds}
        """


checkpoint tribble_reassembly:
    input:
        reads1="<results>/23_tribbles_reassembly/{sample}/trib_reads_R1.fastq.gz",
        reads2="<results>/23_tribbles_reassembly/{sample}/trib_reads_R2.fastq.gz",
        seeds="<results>/23_tribbles_reassembly/{sample}/tribbles_seeds.fasta",
        kmer_params="<results>/20_params_and_stats/200_kmer_params/{sample}.json",
        # kmer_params=outdir / "logs/saute/kmer_params/{sample}.json",
    output:
        # all_vars="<results>/23_tribbles_reassembly/{sample}/all_vars.fasta",
        target_vars="<results>/23_tribbles_reassembly/{sample}/target_vars.fasta",
        graph="<results>/23_tribbles_reassembly/{sample}/graph.gfa",
    params:
        extra="--extend_ends --remove_homopolymer_indels ",
        kmer_threshold=lookup(
            "saute/reassembly/secondary_kmer_threshold", within=pipeline
        ),
        max_var=lookup(query="sample=='{sample}'", cols="homologs", within=sample_table),
        target_cov=lookup("saute/reassembly/target_cov", within=pipeline),
        k1rescale=lookup("saute/reassembly/k1_rescaling", within=pipeline),
        kmers=saute_kmer_expl,
    log:
        "<logs>/<stage>/23_tribbles_reassembly/233_saute/{sample}.log",
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
        --selected_variants {output.target_vars} > {log} 2>&1 || \
        touch {output.graph} {output.target_vars}
        """


rule normal_vars_check:
    input:
        "<results>/22_homologs_assembly/{sample}/tribble_vars.fasta",
    output:
        touch("<results>/23_tribbles_reassembly/{sample}/all_normal.chkp"),


# All probes are normal, no need to do reassembly.
def all_normal(wildcards):
    out_trib = checkpoints.split_saute_assembly.get(sample=wildcards.sample).output.tribbles
    return Path(out_trib).stat().st_size == 0


# Reassembly yield nothing, so take back the original results.
def empty_tribble_asm(wildcards):
    out_trib = checkpoints.tribble_reassembly.get(
        sample=wildcards.sample
    ).output.target_vars
    return Path(out_trib).stat().st_size == 0


rule collect_saute_assemblies:
    input:
        normal="<results>/22_homologs_assembly/{sample}/normal_vars.fasta",
        expl=branch(
            all_normal,
            then="<results>/23_tribbles_reassembly/{sample}/all_normal.chkp",
            otherwise=branch(
                not reasm or empty_tribble_asm,
                then="<results>/22_homologs_assembly/{sample}/tribble_vars.fasta",
                otherwise="<results>/23_tribbles_reassembly/{sample}/target_vars.fasta",
            ),
        ),
    output:
        "<results>/24_postprocessing/240_merging/{sample}.fasta",
    pathvars:
        results="<workdir>/2_saute",
    shell:
        "cat {input.normal} {input.expl} > {output}"


rule cap_tribble_variants:
    input:
        "<results>/24_postprocessing/240_merging/{sample}.fasta",
    output:
        tribble_stats="<results>/24_postprocessing/241_capping/2410_tribbles/{sample}.tsv",
        normal="<results>/24_postprocessing/241_capping/2411_sequences/{sample}.fasta",
    params:
        mode="cap",
        max_var=lookup(query="sample=='{sample}'", cols="homologs", within=sample_table),
        pattern=TARGET_COLLAPSE_PAT,
    log:
        "<logs>/<stage>/24_postprocessing/241_capping/{sample}.log",
    script:
        "../../../src/var_asm_parser.py"


rule fix_homologs_header:
    input:
        "<results>/24_postprocessing/241_capping/2411_sequences/{sample}.fasta",
    output:
        "<results>/24_postprocessing/242_reheading/{sample}.fasta",
    params:
        pattern=SAUTE_PRE_FIX_PAT,
        sample=lambda wildcards: wildcards.sample,
    script:
        "../../../src/fix_headers.py"
