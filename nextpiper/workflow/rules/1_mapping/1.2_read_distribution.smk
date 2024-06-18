targets.append(expand(outdir / "mapped/per_probe/{samples}", samples=sample_list))


# The checkpoint Keyword indicates snakemake to recompute the DAG at this stage.
# This allows the handling of files that are only known at run time.
checkpoint distribute_reads:
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
