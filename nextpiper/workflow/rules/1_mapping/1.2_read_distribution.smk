from reads_distributor import distribute_reads

targets.append(expand(outdir / "mapped_per_probe/{samples}", samples=sample_list))


rule distribute_reads:
    input:
        inbam=outdir / "mapped/{sample}.bam",
    output:
        out=directory(outdir / "mapped_per_probe/{sample}"),
    log:
        outdir / "logs/distribution/{sample}.log",
    threads: 1
    run:
        distribute_reads(input.inbam, output.out)
