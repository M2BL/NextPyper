rule gather_matching_probes:
    input:
        probes=outdir / "translated_probes/longest_cds.fasta",
        tables=expand(
            outdir / "assembled/filtering/seeds_filt_tables/{sample}.tsv",
            sample=sample_list,
        ),
    output:
        outdir / "homolog_prospection/matching_probes.fasta",
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """
        cat {input.tables} | cut -f 2 | sort | uniq > probe_ids.txt
        seqkit grep -nf probe_ids.txt {input.probes} > {output}
        rm probe_ids.txt
        """


use rule make_mmseqs_probes_db as make_mmseqs_matching_probes_db with:
    input:
        outdir / "homolog_prospection/matching_probes.fasta",
    output:
        outdir / "homolog_prospection/homologs_filtering/dbs/matching_probes",
    log:
        outdir / "logs/homolog_prospection/homologs_filtering/make_probes_db.log",


use rule seeds_to_probes_matching as homologs_to_probes_matching with:
    input:
        probes=outdir / "homolog_prospection/homologs_filtering/dbs/matching_probes",
        query=outdir / "saute/target_assembly/{sample}/fixed_vars.fasta",
    output:
        outdir / "homolog_prospection/homologs_filtering/matching_tables/{sample}.tsv",
    log:
        outdir / "logs/homolog_prospection/homologs_filtering/mmseqs/{sample}.log",


use rule seeds_coverage as homologs_coverage with:
    input:
        scfs=outdir / "saute/target_assembly/{sample}/fixed_vars.fasta",
        clean1=outdir / "preprocessed/cleaned/{sample}_R1.fastq.gz",
        clean2=outdir / "preprocessed/cleaned/{sample}_R2.fastq.gz",
    output:
        counts=outdir
        / "homolog_prospection/homologs_filtering/coverage/{sample}.counts",
        metabat=outdir
        / "homolog_prospection/homologs_filtering/coverage/{sample}.metabat",
        hist=outdir / "homolog_prospection/homologs_filtering/coverage/{sample}.hist",
    log:
        outdir / "logs/homolog_prospection/homologs_filtering/coverage/{sample}.log",


use rule seeds_filtering as homologs_filtering with:
    input:
        scfs=outdir / "saute/target_assembly/{sample}/fixed_vars.fasta",
        hits=outdir
        / "homolog_prospection/homologs_filtering/matching_tables/{sample}.tsv",
        covs=outdir / "homolog_prospection/homologs_filtering/coverage/{sample}.metabat",
    output:
        scfs=outdir
        / "homolog_prospection/homologs_filtering/filtered_scfs/{sample}.fasta",
        metrics=outdir
        / "homolog_prospection/homologs_filtering/homolog_filt_tables/{sample}.tsv",
    log:
        outdir
        / "logs/homolog_prospection/homologs_filtering/scfs_filtering/{sample}.log",
    params:
        min_cov=lookup("scf_min_cov", within=homologs_filt_params),
        min_idt=lookup("scf_min_idt", within=homologs_filt_params),
        max_gc=lookup("max_gc", within=homologs_filt_params),
        min_gc=lookup("min_gc", within=homologs_filt_params),
        cov_threshold=lookup("cov_threshold", within=homologs_filt_params),
        cov_dynamic_filt=lookup("cov_dynamic_filt", within=homologs_filt_params),
        separate_probes=lambda wildcards: False,
        tag_scfs=lambda wildcards: False,
        qpat=lambda wildcards: SAUTE_POST_FIX_PAT,
        tpat=lambda wildcards: pattern,


rule estimate_divergence:
    input:
        expand(
            outdir
            / "homolog_prospection/homologs_filtering/homolog_filt_tables/{sample}.tsv",
            sample=sample_list,
        ),
    output:
        outdir / "homolog_prospection/region_separation/divergence_thresholds.json",
    log:
        outdir / "logs/homolog_prospection/region_separation/divergence_estimates.tsv",
    params:
        min_idt=lookup("scf_min_idt", within=homologs_filt_params),
        min_cov=div_est_min_cov,
        flattening_prop=div_est_flat_prop,
    script:
        "../../../src/divergence_estimation.py"
