rule make_mmseqs_long_probes_db:
    input:
        outdir / "translated_probes/longest_cds.fasta",
    output:
        outdir / "assembled/filtering/dbs/probes/probes",
    log:
        outdir / "logs/assembled/filtering/dbs/probes.log",
    conda:
        "../../envs/mmseqs2.yaml"
    shell:
        "mmseqs createdb --dbtype 1 {input} {output} > {log} 2>&1 "


use rule make_mmseqs_raw_assembly_dbs as make_mmseqs_homologs_sample_dbs with:
    input:
        outdir / "assembled/prefixed/{sample}.fasta",
    output:
        outdir / "assembled/filtering/dbs/samples/{sample}",
    log:
        outdir / "logs/assembled/filtering/make_sample_db/{sample}.log",


use rule raw_assembly_to_probes_matching as homologs_to_probes_matching with:
    input:
        probes=outdir / "assembled/filtering/dbs/probes/probes",
        query=outdir / "assembled/filtering/dbs/samples/{sample}",
    output:
        outdir / "assembled/filtering/matching_tables/{sample}.tsv",
    params:
        fields=mmseq_fields,
        evalue=mmseq_evalue,
        min_orf_len=min_orf_len,
        sensitivity=mmseq_sens,
    log:
        outdir / "logs/assembled/filtering/mmseqs/{sample}.log",


checkpoint homologs_filtering:
    input:
        scfs=outdir / "assembled/prefixed/{sample}.fasta",
        table=outdir / "assembled/filtering/matching_tables/{sample}.tsv",
    output:
        directory(outdir / "assembled/filtering/filtered_scfs/{sample}"),
    log:
        outdir / "logs/assembled/filtering/scfs_filtering/{sample}.log",
    params:
        min_cov=seeds_scf_min_cov,
        min_idt=seeds_scf_min_idt,
        separate_probes=lambda wildcards: True,
        qpat=lambda wildcards: False,
        tpat=lambda wildcards: pattern,
    script:
        "../../../src/homolog_filtering.py"


## See Rule 3.1 for further explanation
def aggregate_split(wildcards):
    direcs = []
    for sample in sample_list:
        checkpoint_output = checkpoints.homologs_filtering.get(sample=sample).output[0]
        direcs.append(checkpoint_output)
    return direcs


checkpoint done_assembly:
    input:
        aggregate_split,
    output:
        done=touch(outdir / "logs/dones/splitting.done"),
