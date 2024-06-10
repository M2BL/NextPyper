targets.append(expand(outdir / "mapped/per_probe/{samples}", samples=sample_list))


rule distribute_reads:
    input:
        outdir / "mapped/total/{sample}.bam",
    output:
        directory(outdir / "mapped/per_probe/{sample}"),
    log:
        outdir / "logs/mapped/distribution/{sample}.log",
    threads: 1
    conda:
        "../../envs/distribute_reads.yaml"
    script:
        "../../../src/reads_distributor.py"
