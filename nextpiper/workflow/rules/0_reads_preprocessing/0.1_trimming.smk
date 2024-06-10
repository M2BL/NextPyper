targets.append(
    expand(outdir / "trimmed/{samples}_R{dir}.fastq", samples=sample_list, dir=[1, 2])
)
targets.append(expand(outdir / "report/trimmed/{samples}.json", samples=sample_list))


def get_raw_input_fastq_r1(wildcards):
    return sample_dict[wildcards.sample]["path_forward"]


def get_raw_input_fastq_r2(wildcards):
    return sample_dict[wildcards.sample]["path_reverse"]


rule fastp_pe:
    input:
        in1=get_raw_input_fastq_r1,
        in2=get_raw_input_fastq_r2,
    output:
        trim1=outdir / "trimmed/{sample}_R1.fastq",
        trim2=outdir / "trimmed/{sample}_R2.fastq",
        html=outdir / "report/trimmed/{sample}.html",
        json=outdir / "report/trimmed/{sample}.json",
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
