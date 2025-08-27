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


rule split_saute_assembly:
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
        scfs=outdir / "saute/target_assembly/{sample}/expl_vars.fasta"
        reads1=outdir / "preprocessed/cleaned/{sample}_R1.fastq.gz",
        reads2=outdir / "preprocessed/cleaned/{sample}_R2.fastq.gz",
    output:
        reads1=outdir / "saute/expl_assembly/{sample}/expl_R1.fastq.gz",
        reads2=outdir / "saute/expl_assembly/{sample}/expl_R2.fastq.gz",
    params:
        extra="--proper-pairs-only --exclude-supplementary",
    log:
        outdir / "logs/saute/reassembly/read_collection/{sample}.log",
    threads: 4
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """minimap2 -t {threads} -ax sr {input.scfs} {input.reads1} {input.reads2} 2> {log} | \
        samtools sort -n -u -@ {threads} 2>> {log} | samtools fastq -1 {output.reads1} -2 {output.reads2} - 2>> {log}        
        """ 

rule collect_explosive_seeds:
    input:
        seeds=outdir / "saute/seeds/{sample}.fasta",
        expl=outdir / "saute/target_assembly/{sample}/expl_vars.fasta",
    output:
        seeds=outdir / "saute/expl_assembly/{sample}/expl_seeds.fasta",
    log:
        outdir / "logs/saute/reassembly/seed_collection/{sample}.log",
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """
        awk '/>/{match($0,/-(.*?)_EDGE/, m); print "-"m[1]"_"}' {input.expl} | sort | uniq > seqids_{wildcards.sample}.txt
        seqkit grep -rf seqids.txt  {input.seeds} > {output.seeds}
        rm seqids_{wildcards.sample}.txt
        """
        
use rule saute_assembly as explosive_reassembly with:
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
        kmer_threshold=lookup(
            "saute/reassembly/secondary_kmer_threshold", within=pipeline
        ),
        max_var=lookup("saute/reassembly/max_variants", within=pipeline),
        target_cov=lookup("saute/reassembly/target_cov", within=pipeline),
        kmers=saute_kmer,
    log:
        outdir / "logs/saute/reassembly/assembly/{sample}.log",
    threads: 4
    conda:
        "../../envs/saute.yaml"


use rule split_saute_assembly as collapse_alleles_explosive with:
    input:
        outdir / "saute/expl_assembly/{sample}/target_vars.fasta",
    output:
        normal=outdir / "saute/expl_assembly/{sample}/collapsed_vars.fasta",    
    log:
        outdir / "logs/saute/reassembly/split/{sample}.log",
    

## ToDo: Determine if reassembly is optional or not, and make the rule optional accordingly.

rule collect_saute_assemblies:
    input:
        primary=outdir / "saute/target_assembly/{sample}/collapsed_vars.fasta",
        expl=outdir / "saute/expl_assembly/{sample}/collapsed_vars.fasta",   
    output:
        temp(outdir / "saute/merged/{sample}_merged.fasta"),
    shell:
        "cat {input.primary} {input.expl} > {output}"



rule fix_homologs_header:
    input:
        # outdir / "saute/target_assembly/{sample}/collapsed_vars.fasta",
        outdir / "saute/merged/{sample}_merged.fasta",
    output:
        outdir / "saute/merged/{sample}.fasta",
    params:
        pattern=SAUTE_PRE_FIX_PAT,
        sample=lambda wildcards: wildcards.sample,
    script:
        "../../../src/fix_headers.py"
