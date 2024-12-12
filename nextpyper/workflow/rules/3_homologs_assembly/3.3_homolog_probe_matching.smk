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


rule make_mmseqs_homologs_sample_dbs:
    input:
        outdir / "assembled/prefixed/{sample}.fasta",
    output:
        outdir / "assembled/filtering/dbs/samples/{sample}",
    log:
        outdir / "logs/assembled/filtering/make_sample_db/{sample}.log",
    conda:
        "../../envs/mmseqs2.yaml"
    shell:
        "mmseqs createdb --dbtype 2 {input} {output} > {log} 2>&1"


rule homologs_to_probes_matching:
    input:
        probes=outdir / "assembled/filtering/dbs/probes/probes",
        query=outdir / "assembled/filtering/dbs/samples/{sample}",
    output:
        outdir / "assembled/filtering/matching_tables/{sample}.tsv",
    params:
        fields="query,evalue,qstart,qend,qlen,tstart,tend,tlen,theader,gapopen,nident,mismatch",
        evalue="1.000E-06",
        min_orf_len=15,
        sensitivity=7.5,
    log:
        outdir / "logs/assembled/filtering/mmseqs/{sample}.log",
    threads: 4
    conda:
        "../../envs/mmseqs2.yaml"
    shell:
        """
        mkdir -p temp_{wildcards.sample}
        mmseqs search {input.query} {input.probes} {wildcards.sample}_results temp_{wildcards.sample} --threads {threads} -s {params.sensitivity} -e {params.evalue} --min-length {params.min_orf_len} --remove-tmp-files -a > {log} 2>&1
        mmseqs convertalis {input.query} {input.probes} {wildcards.sample}_results {output} --format-mode 4 --format-output {params.fields} --threads {threads} >> {log} 2>&1
        rm -r temp_{wildcards.sample}
        rm *_results.*
        """


checkpoint homologs_filtering:
    input:
        scfs=outdir / "assembled/prefixed/{sample}.fasta",
        table=outdir / "assembled/filtering/matching_tables/{sample}.tsv",
    output:
        directory(outdir / "assembled/filtering/filtered_scfs/{sample}"),
    log:
        outdir / "logs/assembled/filtering/scfs_filtering/{sample}.log",
    params:
        min_cov=homolog_scf_min_cov,
        min_idt=homolog_scf_min_idt,
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
