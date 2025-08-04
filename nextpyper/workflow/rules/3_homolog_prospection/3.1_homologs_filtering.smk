rule gather_matching_probes:
    input:
        probes=outdir / "translated_probes/longest_cds.fasta",
        tables=expand(
            outdir / "logs/assembled/filtering/scfs_filtering/{samples}.log",
            samples=sample_list,
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
        outdir / "homolog_prospection/homologs_filtering/dbs/probes/matching_probes",
    log:
        outdir / "logs/homolog_prospection/homologs_filtering/make_probes_db.log",


use rule seeds_to_probes_matching as homologs_to_probes_matching with:
    input:
        probes=outdir
        / "homolog_prospection/homologs_filtering/dbs/probes/matching_probes",
        query=outdir / "saute/target_assembly/{sample}/fixed_vars.fasta",
    output:
        outdir / "homolog_prospection/homologs_filtering/matching_tables/{sample}.tsv",
    log:
        outdir / "logs/homolog_prospection/homologs_filtering/mmseqs/{sample}.log",


use rule seeds_filtering as homologs_filtering with:
    input:
        scfs=outdir / "saute/target_assembly/{sample}/fixed_vars.fasta",
        table=outdir
        / "homolog_prospection/homologs_filtering/matching_tables/{sample}.tsv",
    output:
        outdir / "homolog_prospection/homologs_filtering/filtered_scfs/{sample}.fasta",
    log:
        outdir
        / "logs/homolog_prospection/homologs_filtering/scfs_filtering/{sample}.log",
    params:
        min_cov=homolog_scf_min_cov,
        min_idt=homolog_scf_min_idt,
        separate_probes=lambda wildcards: False,
        separate_scfs=lambda wildcards: False,
        qpat=lambda wildcards: SAUTE_POST_FIX_PAT,
        tpat=lambda wildcards: pattern,


rule estimate_divergence:
    input:
        expand(
            outdir
            / "logs/homolog_prospection/homologs_filtering/scfs_filtering/{samples}.log",
            samples=sample_list,
        ),
    output:
        outdir / "homolog_prospection/region_separation/divergence_thresholds.json",
    log:
        outdir / "logs/homolog_prospection/region_separation/divergence_estimates.tsv",
    params:
        min_idt=homolog_scf_min_idt,
        min_cov=div_est_min_cov,
        flattening_prop=div_est_flat_prop,
    script:
        "../../../src/divergence_estimation.py"
