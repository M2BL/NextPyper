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
        min_cov=lookup("divergence_estimation/min_cov", within=pipeline),
        flattening_prop=lookup("divergence_estimation/flatenning_prop", within=pipeline),
    script:
        "../../../src/divergence_estimation.py"


rule estimate_intron_ceiling:
    input:
        expand(outdir / "logs/preprocessing/fastp/{sample}.json", sample=sample_list),
    output:
        outdir / "homolog_prospection/region_separation/intron_ceilings.json",
    run:
        max_intron = {
            file.stem: read_insert_size(file) * 2 for file in map(Path, input)
        }
        Path(output[0]).write_text(json.dumps(max_intron, indent=4))


use rule distribute_seeds as per_probe_scaffold_grouping with:
    input:
        expand(
            outdir
            / "homolog_prospection/homologs_filtering/filtered_scfs/{sample}.fasta",
            sample=sample_list,
        ),
    output:
        expand(
            outdir / "homolog_prospection/region_separation/input_scfs/{probe}.fasta",
            probe=probes_list,
        ),
    log:
        outdir / "logs/homolog_prospection/region_separation/scfs_grouping.log",
    params:
        pattern=lambda wildcards: SAUTE_POST_FIX_PAT,


use rule distribute_seeds as split_matching_probes with:
    input:
        probes=outdir / "homolog_prospection/matching_probes.fasta",
        tables=expand(
            outdir
            / "homolog_prospection/homologs_filtering/homolog_filt_tables/{sample}.tsv",
            sample=sample_list,
        ),
    output:
        expand(
            outdir
            / "homolog_prospection/region_separation/input_probes/{probe}.fasta",
            probe=probes_list,
        ),
    log:
        outdir / "logs/homolog_prospection/region_separation/probe_grouping.log",
    params:
        pattern=lambda wildcards: probe_pattern,
        mode="multi_probes" if multi_probes else "single_probes",


rule separate_cds_by_regions:
    input:
        probes=outdir
        / "homolog_prospection/region_separation/input_probes/{probe}.fasta",
        scfs=outdir / "homolog_prospection/region_separation/input_scfs/{probe}.fasta",
        div_map=outdir
        / "homolog_prospection/region_separation/divergence_thresholds.json",
        max_intron_map=outdir
        / "homolog_prospection/region_separation/intron_ceilings.json",
    output:
        directory(
            outdir
            / "homolog_prospection/region_separation/separation_output/scfs/{probe}"
        ),
    params:
        force_global_idt=lookup("enforce_global_idt_threshold", within=reg_sep),
        min_global_identity=lookup("min_global_identity", within=reg_sep),
        min_fragment_cov=lookup("min_fragment_cov", within=reg_sep),
        min_exonic_length=lookup("min_exonic_length", within=reg_sep),
        max_intron_length=lookup("max_intron_length", within=reg_sep),
        substitution_matrix=blosum62,
    log:
        outdir / "logs/homolog_prospection/region_separation/separation/{probe}.log",
    threads: 2
    conda:
        "../../envs/clustering.yaml"
    script:
        "../../../src/miniprot.py"


rule align_regions:
    input:
        outdir / "homolog_prospection/region_separation/separation_output/scfs/{probe}",
    output:
        directory(
            outdir / "homolog_prospection/region_separation/alns/{probe}",
        ),
    params:
        lookup("mafft_params", within=pipeline),
    log:
        outdir / "logs/homolog_prospection/region_separation/alns/{probe}.log",
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
            scfs=expand(
                outdir
                / "homolog_prospection/region_separation/separation_output/scfs/{probe}",
                probe=probes_list,
            ),
            tribbles=expand(
                outdir / "logs/saute/tribbles/{sample}.tsv", sample=sample_list
            ),
            ## Disabled temporarily. See chimera_tagging rule in 3.2 for details.
            # chimera_tags=expand(
            #     outdir
            #     / "homolog_prospection/homologs_filtering/chimera_tagging/{sample}.tsv",
            #     sample=sample_list,
            # ),
        output:
            expand(
                outdir
                / f"homolog_prospection/region_separation/consolidated/per_sample/{kind}/{{sample}}.fasta",
                sample=sample_list,
            ),
        log:
            outdir
            / f"logs/homolog_prospection/region_separation/consolidated/grouping/{kind}.log",
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
            scfs=outdir
            / f"homolog_prospection/region_separation/consolidated/per_sample/{kind}/{{sample}}.fasta",
            clean1=outdir / "preprocessed/cleaned/{sample}_R1.fastq.gz",
            clean2=outdir / "preprocessed/cleaned/{sample}_R2.fastq.gz",
        output:
            counts=outdir
            / f"homolog_prospection/region_separation/consolidated/coverage/{kind}/{{sample}}.counts",
            metabat=outdir
            / f"homolog_prospection/region_separation/consolidated/coverage/{kind}/{{sample}}.metabat",
        params:
            extra="--proper-pairs-only --exclude-supplementary",
        log:
            outdir
            / f"logs/homolog_prospection/region_separation/consolidated/coverage/{kind}/{{sample}}.log",
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
