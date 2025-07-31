checkpoint per_probe_scaffold_grouping:
    input:
        expand(
            outdir / "homolog_prospection/allele_collapsing/{sample}.fasta",
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
        probes=probes_list,
        mode="scfs",
    script:
        "../../../src/multi_seq_probes.py"


checkpoint split_matching_probes:
    input:
        probes=outdir / "homolog_prospection/matching_probes.fasta",
        tables=expand(
            outdir
            / "logs/homolog_prospection/homologs_filtering/scfs_filtering/{samples}.log",
            samples=sample_list,
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
        pattern=lambda wildcards: pattern,
        probes=probes_list,
        mode="multi_probes" if multi_probes else "single_probes",
    conda:
        "../../envs/preprocessing.yaml"
    script:
        "../../../src/multi_seq_probes.py"


checkpoint separate_cds_by_regions:
    input:
        probes=outdir
        / "homolog_prospection/region_separation/input_probes/{probe}.fasta",
        scfs=outdir / "homolog_prospection/region_separation/input_scfs/{probe}.fasta",
        div_map=outdir
        / "homolog_prospection/region_separation/divergence_thresholds.json",
    output:
        directory(
            outdir
            / "homolog_prospection/region_separation/separation_output/scfs/{probe}"
        ),
    params:
        min_fragment_cov=min_fragment_cov,
        min_exonic_length=min_exonic_length,
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
        "--auto --reorder",
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
