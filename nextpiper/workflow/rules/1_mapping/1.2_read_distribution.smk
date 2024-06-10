from reads_distributor import distribute_reads

targets.append(expand(outdir / "mapped/per_probe/{samples}", samples=sample_list))


rule distribute_reads:
    input:
        inbam=outdir / "mapped/total/{sample}.bam",
    output:
        out=directory(outdir / "mapped/per_probe/{sample}"),
    log:
        outdir / "logs/mapped/distribution/{sample}.log",
    threads: 1
    run:
        distribute_reads(input.inbam, output.out)
