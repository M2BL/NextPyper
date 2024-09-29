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
        trim1=outdir / "preprocessed/trimmed/{sample}_R1.fastq",
        trim2=outdir / "preprocessed/trimmed/{sample}_R2.fastq",
        html=outdir / "logs/trimmed/{sample}.html",
        json=outdir / "logs/trimmed/{sample}.json",
    log:
        outdir / "logs/fastp/{sample}.log",
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
        in1=outdir / "preprocessed/trimmed/{sample}_R1.fastq",
        in2=outdir / "preprocessed/trimmed/{sample}_R2.fastq",
        ref=probes_path.resolve(),
    output:
        outm1=outdir / "preprocessed/filtered/{sample}_R1.fastq",
        outm2=outdir / "preprocessed/filtered/{sample}_R2.fastq",
    log:
        outdir / "logs/preprocessing/bbduk/{sample}.log",
    params:
        others="-Xmx4g",
        command="bbduk.sh",
        k=19,
    threads: 4
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        "(bbduk.sh "
        "{params.others} "
        "in1={input.in1} "
        "in2={input.in2} "
        "ref={input.ref} "
        "outm1={output.outm1} "
        "outm2={output.outm1} "
        "k={params.k} "
        "threads={threads} ) >> {log} 2>&"


# rule matching_probes:
#     input:
#         input=[
#             outdir / "preprocessed/trimmed/{sample}_R1.fastq",
#             outdir / "preprocessed/trimmed/{sample}_R2.fastq",
#         ],
#         ref=probes_path.resolve(),
#     output:
#         outm=[
#             outdir / "preprocessed/filtered/{sample}_R1.fastq",
#             outdir / "preprocessed/filtered/{sample}_R2.fastq",
#         ],
#     log:
#         outdir / "logs/preprocessing/bbduk/{sample}.log",
#     params:
#         Xmx4g,
#         command="bbduk.sh",
#         k=19,
#     threads: 4
#     wrapper:
#         "v4.3.0/bio/bbtools"
