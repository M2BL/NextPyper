# The input function map the sample at hand to its input files (specified in the samples table):
def get_raw_input_fastq_r1(wildcards):
    return sample_dict[wildcards.sample]["path_forward"]


def get_raw_input_fastq_r2(wildcards):
    return sample_dict[wildcards.sample]["path_reverse"]


rule fastp_pe:
    input:
        in1=get_raw_input_fastq_r1,
        in2=get_raw_input_fastq_r2,
    output:
        trim1=outdir / "preprocessed/trimmed/{sample}_R1.fastq.gz",
        trim2=outdir / "preprocessed/trimmed/{sample}_R2.fastq.gz",
        html=outdir / "logs/preprocessing/fastp/{sample}.html",
        json=outdir / "logs/preprocessing/fastp/{sample}.json",
    log:
        outdir / "logs/preprocessing/fastp/{sample}.log",
    params:
        extra="--trim_poly_g --trim_poly_x --low_complexity_filter --cut_right",
    threads: 4
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        "(fastp --thread {threads} "
        "{params.extra} "
        "--in1 {input.in1} --in2 {input.in2} "
        "--out1 {output.trim1} --out2 {output.trim2} "
        "--html {output.html} "
        "--json {output.json} ) 2> {log} "


rule matching_probes:
    input:
        in1=outdir / "preprocessed/trimmed/{sample}_R1.fastq.gz",
        in2=outdir / "preprocessed/trimmed/{sample}_R2.fastq.gz",
        ref=silva_db,
    output:
        out1=outdir / "preprocessed/cleaned/{sample}_R1.fastq",
        out2=outdir / "preprocessed/cleaned/{sample}_R2.fastq",
    log:
        outdir / "logs/preprocessing/bbduk_cleaning/{sample}.log",
    params:
        others=other_bbduk,
        k=19,
    threads: 4
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        "(bbduk.sh {params.others} "
        "in1={input.in1} "
        "in2={input.in2} "
        "ref={input.ref} "
        "out1={output.out1} "
        "out2={output.out2} "
        "k={params.k} "
        "threads={threads} ) 2> {log} "


rule matching_probes:
    input:
        in1=outdir / "preprocessed/cleaned/{sample}_R1.fastq.gz",
        in2=outdir / "preprocessed/cleaned/{sample}_R2.fastq.gz",
        ref=probes_path.resolve(),
    output:
        outm1=outdir / "preprocessed/filtered/{sample}_R1.fastq",
        outm2=outdir / "preprocessed/filtered/{sample}_R2.fastq",
    log:
        outdir / "logs/preprocessing/bbduk/{sample}.log",
    params:
        others=other_bbduk,
        k=bbduk_k,
    threads: 4
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        "(bbduk.sh {params.others} "
        "in1={input.in1} "
        "in2={input.in2} "
        "ref={input.ref} "
        "outm1={output.outm1} "
        "outm2={output.outm2} "
        "k={params.k} "
        "threads={threads} ) 2> {log} "
