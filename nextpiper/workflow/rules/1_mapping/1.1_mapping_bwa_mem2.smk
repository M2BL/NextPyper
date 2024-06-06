targets.append(expand(outdir / "mapped/{samples}.bam", samples=sample_list))


rule probes_symlink:
    input:
        probes.resolve(),
    output:
        outdir / "mapped/map_index/probes.fasta",
    shell:
        "ln -s {input} {output}"


rule bwa_mem2_index:
    input:
        rules.probes_symlink.output,
    output:
        multiext(
            str(outdir / "mapped/map_index/probes.fasta"),
            ".0123",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
        ),
    log:
        outdir / "logs/bwa-mem2_index/probes.log",
    wrapper:
        "v3.11.0/bio/bwa-mem2/index"


rule bwa_mem2_mem:
    input:
        reads=[
            outdir / "trimmed/{sample}_R1.fastq",
            outdir / "trimmed/{sample}_R2.fastq",
        ],
        idx=rules.bwa_mem2_index.output,
    output:
        outdir / "mapped/{sample}.bam",
    log:
        outdir / "logs/bwa_mem2/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",
        sort_order="queryname",
        sort_extra="",
    threads: 8
    wrapper:
        "v3.11.0/bio/bwa-mem2/mem"
