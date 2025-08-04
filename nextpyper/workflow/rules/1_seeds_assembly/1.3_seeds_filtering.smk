use rule raw_assembly_to_probes_matching as seeds_to_probes_matching with:
    input:
        probes=outdir / "assembled/filtering/dbs/probes/probes",
        query=outdir / "assembled/prefixed/{sample}.fasta",
    output:
        outdir / "assembled/filtering/matching_tables/{sample}.tsv",
    params:
        fields=mmseq_fields,
        evalue=mmseq_evalue,
        min_orf_len=min_orf_len,
        sensitivity=mmseq_sens,
    log:
        outdir / "logs/assembled/filtering/mmseqs/{sample}.log",


rule seeds_filtering:
    input:
        scfs=outdir / "assembled/prefixed/{sample}.fasta",
        table=outdir / "assembled/filtering/matching_tables/{sample}.tsv",
    output:
        outdir / "assembled/filtering/filtered_scfs/{sample}.fasta",
    log:
        outdir / "logs/assembled/filtering/scfs_filtering/{sample}.log",
    params:
        min_cov=seeds_scf_min_cov,
        min_idt=seeds_scf_min_idt,
        separate_probes=lambda wildcards: False,
        tag_scfs=lambda wildcards: True,
        qpat=lambda wildcards: False,
        tpat=lambda wildcards: pattern,
    script:
        "../../../src/homolog_filtering.py"
