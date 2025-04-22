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
        trim_qual=trim_qual,
        trim_min_len=trim_min_len,
        extra="--trim_poly_g --trim_poly_x --low_complexity_filter --cut_tail",
    threads: 4
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        "(fastp --thread {threads} "
        "{params.extra} "
        "--cut_mean_quality {params.trim_qual} "
        "--length_required {params.trim_min_len} "
        "--in1 {input.in1} --in2 {input.in2} "
        "--out1 {output.trim1} --out2 {output.trim2} "
        "--html {output.html} "
        "--json {output.json} ) 2> {log} "
