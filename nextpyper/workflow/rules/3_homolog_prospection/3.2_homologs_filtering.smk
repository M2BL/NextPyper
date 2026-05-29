pathvars:
    stage="3_homolog_prospection/31_homologs_filtering",
    results="<workdir>/<stage>",
    inter="3_homolog_prospection/30_allele_collapsing",


rule gather_matching_probes:
    input:
        probes="<workdir>/0_preprocessing/longest_cds.fasta",
        tables=expand(
            "<workdir>/1_assembly/12_filtering/122_filt_seeds/1220_summary_tables/{sample}.tsv",
            sample=sample_list,
        ),
    output:
        "<workdir>/3_homolog_prospection/matching_probes.fasta",
    shadow:
        "shallow"
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """
        cat {input.tables} | cut -f 2 | sort | uniq > probe_ids.txt
        seqkit grep -f probe_ids.txt {input.probes} > {output}
        """


use rule make_diamond_probes_db as make_diamond_matching_probes_db with:
    input:
        "<workdir>/3_homolog_prospection/matching_probes.fasta",
    output:
        "<results>/matching_probes.dmnd",
    log:
        "<logs>/<stage>/make_probes_db.log",


use rule raw_assembly_to_probes_matching as homologs_to_probes_matching with:
    input:
        probes="<results>/matching_probes.dmnd",
        query="<workdir>/<inter>/{sample}.fasta",
    output:
        "<results>/310_matching_tables/{sample}.tsv",
    log:
        "<logs>/<stage>/310_matching_diamond/{sample}.log",


use rule seeds_coverage as homologs_coverage with:
    input:
        scfs="<workdir>/<inter>/{sample}.fasta",
        clean1="<workdir>/0_preprocessing/02_cleaning/{sample}_R1.fastq.gz",
        clean2="<workdir>/0_preprocessing/02_cleaning/{sample}_R2.fastq.gz",
    output:
        counts="<results>/311_coverage/{sample}.counts",
        metabat="<results>/311_coverage/{sample}.metabat",
    log:
        "<logs>/<stage>/311_coverage/{sample}.log",


use rule seeds_filtering as homologs_filtering with:
    input:
        scfs="<workdir>/<inter>/{sample}.fasta",
        hits="<results>/310_matching_tables/{sample}.tsv",
        covs="<results>/311_coverage/{sample}.metabat",
    output:
        metrics="<results>/312_filt_homologs/3120_summary_tables/{sample}.tsv",
        scfs="<results>/312_filt_homologs/3121_sequences/{sample}.fasta",
    log:
        "<logs>/<stage>/312_filt_homologs/{sample}.log",
    params:
        min_cov=lookup("scf_min_cov", within=homologs_filt_params),
        min_idt=lookup("scf_min_idt", within=homologs_filt_params),
        max_gc=lookup("max_gc", within=homologs_filt_params),
        min_gc=lookup("min_gc", within=homologs_filt_params),
        cov_threshold=lookup("cov_threshold", within=homologs_filt_params),
        cov_dynamic_filt=lookup("cov_dynamic_filt", within=homologs_filt_params),
        tag_scfs=lambda wildcards: False,
        qpat=lambda wildcards: SAUTE_POST_FIX_PAT,


## As of v2.30.2, vsearch segfaults unpredictably with some inputs, crashing the whole pipeline.
## This feature is disabled until a fix is available.
## ToDo: A fix is available since v2.30.3. Test and consider re-enabling chimera tagging.
# rule chimera_tagging:
#     input:
#         outdir / "homolog_prospection/homologs_filtering/filtered_scfs/{sample}.fasta",
#     output:
#         outdir / "homolog_prospection/homologs_filtering/chimera_tagging/{sample}.tsv",
#     log:
#         outdir
#         / "logs/homolog_prospection/homologs_filtering/chimera_tagging/{sample}.log",
#     conda:
#         "../../envs/clustering.yaml"
#     shell:
#         """
#         vsearch --chimeras_denovo {input} \
#             --tabbedout {output} 2> {log}
#         """
