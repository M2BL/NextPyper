targets.append(
    expand(outdir / "logs/assembled/collect/{samples}.chkpt", samples=sample_list)
)


rule bam2fastq:
    input:
        aux=rules.distribute_reads.output,
        bam=outdir / "mapped/per_probe/{sample}/{probe}.bam",
    output:
        out1=outdir / "assembled/inputs/{sample}/{probe}_R1.fastq",
        out2=outdir / "assembled/inputs/{sample}/{probe}_R2.fastq",
    log:
        outdir / "logs/assembled/inputs/{sample}/{probe}.log",
    conda:
        "../../envs/assembly_spades.yaml"
    shell:
        "samtools fastq -1 {output.out1} -2 {output.out2} {input.bam} 2> {log}"


rule spades_assembly:
    input:
        in1=outdir / "assembled/inputs/{sample}/{probe}_R1.fastq",
        in2=outdir / "assembled/inputs/{sample}/{probe}_R2.fastq",
    output:
        out_dir=directory(outdir / "assembled/spades/{sample}/{probe}"),
        contigs=outdir / "assembled/spades/{sample}/{probe}/contigs.fasta",
        gfa=outdir
        / "assembled/spades/{sample}/{probe}/assembly_graph_with_scaffolds.gfa",
    params:
        "",
    log:
        outdir / "logs/assembled/spades/{sample}/{probe}.log",
    threads: 1
    conda:
        "../../envs/assembly_spades.yaml"
    shell:
        "spades.py -t {threads} {params} -1 {input.in1} -2 {input.in2} -o {output.out_dir} > {log} 2>&1"


# Input function to handle dynamically generated files via checkpoints.
# See: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
def aggregate_asms(wildcards):
    checkpoint_output = checkpoints.distribute_reads.get(**wildcards).output[0]
    return expand(
        outdir / "assembled/spades/{sample}/{probe}",
        sample=wildcards.sample,
        probe=glob_wildcards(os.path.join(checkpoint_output, "{probe}.bam")).probe,
    )  # glob_wilcards(path, wildcard_pattern) tries to match the given wildcard to
    # the closest file generated so far by the workflow.


# Collect all the dinamycally generated files (at run time)
rule collect_assemblies:
    input:
        aggregate_asms,
    output:
        outdir / "logs/assembled/collect/{sample}.chkpt",
    shell:
        "echo {input} | tr '[:space:]' '\n' >> {output}"
