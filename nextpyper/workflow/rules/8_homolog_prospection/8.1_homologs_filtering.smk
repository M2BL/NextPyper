def gather_tables(wildcards):
    return expand(
        outdir / "assembled/filtering/matching_tables/{samples}.tsv",
        samples=sample_list,
    )


rule gather_matching_probes:
    input:
        probes=outdir / "translated_probes/longest_cds.fasta",
        tables=expand(
            outdir / "assembled/filtering/matching_tables/{samples}.tsv",
            samples=sample_list,
        ),
    output:
        outdir / "homolog_prospection/matching_probes.fasta",
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """
        cat {input.tables} | cut -f 9 | sort | uniq > probe_ids.txt
        seqkit grep -nf probe_ids.txt {input.probes} > {output}
        rm probe_ids.txt
        """


rule make_mmseqs_probe_db:
    input:
        outdir / "homolog_prospection/matching_probes.fasta",
    output:
        outdir / "homolog_prospection/candidates_filtering/dbs/probes/matching_probes",
    log:
        outdir / "logs/homolog_prospection/candidates_filtering/make_probes_db.log",
    conda:
        "../../envs/mmseqs2.yaml"
    shell:
        "mmseqs createdb --dbtype 1 {input} {output} > {log} 2>&1"


use rule make_mmseqs_homologs_sample_dbs as make_mmseqs_sample_dbs with:
    input:
        outdir / "saute/target_assembly/{sample}/target_vars.fasta",
    output:
        outdir / "homolog_prospection/candidates_filtering/dbs/samples/{sample}",
    log:
        outdir
        / "logs/homolog_prospection/candidates_filtering/make_sample_db/{sample}.log",


use rule homologs_to_probes_matching as candidates_to_probes_matching with:
    input:
        probes=outdir
        / "homolog_prospection/candidates_filtering/dbs/probes/matching_probes",
        query=outdir / "homolog_prospection/candidates_filtering/dbs/samples/{sample}",
    output:
        outdir / "homolog_prospection/candidates_filtering/matching_tables/{sample}.tsv",
    log:
        outdir / "logs/homolog_prospection/candidates_filtering/mmseqs/{sample}.log",


use rule homologs_filtering as candidates_filtering with:
    input:
        scfs=outdir / "saute/target_assembly/{sample}/target_vars.fasta",
        table=outdir
        / "homolog_prospection/candidates_filtering/matching_tables/{sample}.tsv",
    output:
        outdir / "homolog_prospection/candidates_filtering/filtered_scfs/{sample}.fasta",
    log:
        outdir
        / "logs/homolog_prospection/candidates_filtering/scfs_filtering/{sample}.log",
    params:
        min_cov=homolog_scf_min_cov,
        min_idt=homolog_scf_min_idt,
        separate_probes=lambda wildcards: False,
        tpat=lambda wildcards: pattern,
