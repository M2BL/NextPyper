use rule raw_assembly_to_probes_matching as seeds_to_probes_matching with:
    input:
        probes=outdir / "assembled/filtering/probes.dmnd",
        query=outdir / "assembled/extension/{sample}.fasta",
    output:
        outdir / "assembled/filtering/matching_tables/{sample}.tsv",
    log:
        outdir / "logs/assembled/filtering/diamond/{sample}.log",


rule seeds_coverage:
    input:
        scfs=outdir / "assembled/extension/{sample}.fasta",
        clean1=outdir / "preprocessed/cleaned/{sample}_R1.fastq.gz",
        clean2=outdir / "preprocessed/cleaned/{sample}_R2.fastq.gz",
    output:
        counts=outdir / "assembled/filtering/coverage/{sample}.counts",
        metabat=outdir / "assembled/filtering/coverage/{sample}.metabat",
    params:
        extra="--proper-pairs-only --exclude-supplementary",
    log:
        outdir / "logs/assembled/filtering/coverage/{sample}.log",
    threads: 4
    shadow:
        "shallow"
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """
        minimap2 -t {threads} -ax sr {input.scfs} {input.clean1} {input.clean2} 2> {log} | \
        samtools sort -u -@ {threads} > tmp_{wildcards.sample}.bam 2>> {log}
        coverm contig -m count {params.extra} -b tmp_{wildcards.sample}.bam  > {output.counts} 2>> {log}
        coverm contig -m metabat {params.extra} -b tmp_{wildcards.sample}.bam  > {output.metabat} 2>> {log}
        """


rule seeds_filtering:
    input:
        scfs=outdir / "assembled/extension/{sample}.fasta",
        hits=outdir / "assembled/filtering/matching_tables/{sample}.tsv",
        covs=outdir / "assembled/filtering/coverage/{sample}.metabat",
    output:
        scfs=outdir / "assembled/filtering/filtered_scfs/{sample}.fasta",
        metrics=outdir / "assembled/filtering/seeds_filt_tables/{sample}.tsv",
    log:
        outdir / "logs/assembled/filtering/seeds_filtering/{sample}.log",
    params:
        min_cov=lookup("scf_min_cov", within=seeds_filt_params),
        min_idt=lookup("scf_min_idt", within=seeds_filt_params),
        max_gc=lookup("max_gc", within=seeds_filt_params),
        min_gc=lookup("min_gc", within=seeds_filt_params),
        cov_threshold=lookup("cov_threshold", within=seeds_filt_params),
        cov_dynamic_filt=lookup("cov_dynamic_filt", within=seeds_filt_params),
        separate_probes=lambda wildcards: False,
        tag_scfs=lambda wildcards: False,
        qpat=lambda wildcards: SEED_PAT,
        tpat=lambda wildcards: probe_pattern,
    script:
        "../../../src/homolog_filtering.py"
