pathvars:
    stage="3_homolog_prospection/32_region_separation",
    results="<workdir>/<stage>",
    inter="3_homolog_prospection/31_homologs_filtering/312_filt_homologs",


rule estimate_divergence:
    input:
        expand("<workdir>/<inter>/3120_summary_tables/{sample}.tsv", sample=sample_list),
    output:
        "<results>/divergence_thresholds.json",
    log:
        "<logs>/<stage>/divergence_estimates.tsv",
    params:
        min_idt=lookup("scf_min_idt", within=homologs_filt_params),
        min_cov=lookup("divergence_estimation/min_cov", within=pipeline),
        flattening_prop=lookup("divergence_estimation/flatenning_prop", within=pipeline),
    script:
        "../../../src/divergence_estimation.py"


rule estimate_intron_ceiling:
    input:
        expand("<logs>/0_preprocessing/01_trimming_fastp/{sample}.json", sample=sample_list),
    output:
        "<results>/intron_ceilings.json",
    run:
        max_intron = {
            file.stem: read_insert_size(file) * 2 for file in map(Path, input)
        }
        Path(output[0]).write_text(json.dumps(max_intron, indent=4))


use rule distribute_seeds as per_probe_scaffold_grouping with:
    input:
        expand("<workdir>/<inter>/3121_sequences/{sample}.fasta", sample=sample_list),
    output:
        expand("<results>/320_inputs/3201_scfs/{probe}.fasta", probe=probes_list),
    log:
        "<logs>/<stage>/scfs_grouping.log",
    params:
        pattern=lambda wildcards: SAUTE_POST_FIX_PAT,


use rule distribute_seeds as split_matching_probes with:
    input:
        probes="<workdir>/3_homolog_prospection/matching_probes.fasta",
        tables=expand("<workdir>/<inter>/3120_summary_tables/{sample}.tsv", sample=sample_list),
    output:
        expand("<results>/320_inputs/3200_probes/{probe}.fasta", probe=probes_list),
    log:
        "<logs>/<stage>/probe_grouping.log",
    params:
        pattern=lambda wildcards: probe_pattern,
        mode="multi_probes" if multi_probes else "single_probes",


rule separate_cds_by_regions:
    input:
        probes="<results>/320_inputs/3200_probes/{probe}.fasta",
        scfs="<results>/320_inputs/3201_scfs/{probe}.fasta",
        div_map="<results>/divergence_thresholds.json",
        max_intron_map="<results>/intron_ceilings.json",
    output:
        directory("<results>/321_output/3211_scfs/{probe}"),
    params:
        force_global_idt=lookup("enforce_global_idt_threshold", within=reg_sep),
        min_global_identity=lookup("min_global_identity", within=reg_sep),
        min_fragment_cov=lookup("min_fragment_cov", within=reg_sep),
        min_exonic_length=lookup("min_exonic_length", within=reg_sep),
        max_intron_length=lookup("max_intron_length", within=reg_sep),
        substitution_matrix=blosum62,
        probes_outdir=outdir / "workflow/3_homolog_prospection/32_region_separation/321_output/3210_probes",
    log:
        "<logs>/<stage>/321_miniprot/{probe}.log",
    threads: 2
    conda:
        "../../envs/clustering.yaml"
    script:
        "../../../src/miniprot.py"


rule align_regions:
    input:
        "<results>/321_output/3211_scfs/{probe}",
    output:
        directory("<workdir>/3_homolog_prospection/33_alns/{probe}"),
    params:
        lookup("mafft_params", within=pipeline),
    log:
        "<logs>/3_homolog_prospection/33_alns/{probe}.log",
    threads: 2
    conda:
        "../../envs/alignment.yaml"
    shell:
        """
        rm -f {log}
        mkdir -p {output}

        for file in $(find {input} -name "*.fasta"); do
            name=$(basename $file)
            nseqs=$(grep -c "^>" $file)

            if [ "$nseqs" -gt 1 ]; then
                mafft --thread {threads} {params} $file > {output}/$name 2>> {log}
            else
                cp $file {output}/$name
            fi 
        done
        """


for kind in ("exons", "genetigs", "supercontigs"):

    rule:
        name:
            f"collect_{kind}"
        input:
            scfs=expand("<results>/321_output/3211_scfs/{probe}", probe=probes_list),
            tribbles=expand("<workdir>/2_saute/24_postprocessing/241_capping/2410_tribbles/{sample}.tsv", sample=sample_list),
            ## Disabled temporarily. See chimera_tagging rule in 3.2 for details.
            # chimera_tags=expand(
            #     outdir
            #     / "homolog_prospection/homologs_filtering/chimera_tagging/{sample}.tsv",
            #     sample=sample_list,
            # ),
        output:
            expand(f"<workdir>/3_homolog_prospection/34_collating_homologs/340_per_sample/{kind}/{{sample}}.fasta", sample=sample_list),
        log:
            f"<logs>/3_homolog_prospection/34_collating_homologs/340_grouping/{kind}.log",
        params:
            pattern=lambda wildcards: COMP_FINAL_PAT,
            probes=probes_list,
            mode=kind,
        script:
            "../../../src/multi_seq_probes.py"

    rule:
        name:
            f"{kind}_coverage"
        input:
            scfs=f"<workdir>/3_homolog_prospection/34_collating_homologs/340_per_sample/{kind}/{{sample}}.fasta",
            clean1="<workdir>/0_preprocessing/02_cleaning/{sample}_R1.fastq.gz",
            clean2="<workdir>/0_preprocessing/02_cleaning/{sample}_R2.fastq.gz",
        output:
            counts=f"<workdir>/3_homolog_prospection/34_collating_homologs/341_coverage/{kind}/{{sample}}.counts",
            metabat=f"<workdir>/3_homolog_prospection/34_collating_homologs/341_coverage/{kind}/{{sample}}.metabat",
        params:
            extra="--proper-pairs-only --exclude-supplementary",
        log:
            f"<logs>/3_homolog_prospection/34_collating_homologs/341_coverage/{kind}/{{sample}}.log",
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
