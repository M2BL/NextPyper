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


rule make_mmseqs_sample_dbs:
    input:
        outdir / "saute/target_assembly/{sample}/target_vars.fasta",
    output:
        outdir / "homolog_prospection/candidates_filtering/dbs/samples/{sample}",
    log:
        outdir
        / "logs/homolog_prospection/candidates_filtering/make_sample_db/{sample}.log",
    conda:
        "../../envs/mmseqs2.yaml"
    shell:
        "mmseqs createdb --dbtype 2 {input} {output} > {log} 2>&1"


rule candidates_to_probes_matching:
    input:
        probes=outdir
        / "homolog_prospection/candidates_filtering/dbs/probes/matching_probes",
        query=outdir / "homolog_prospection/candidates_filtering/dbs/samples/{sample}",
    output:
        outdir / "homolog_prospection/candidates_filtering/matching_tables/{sample}.tsv",
    params:
        fields="query,evalue,qstart,qend,qlen,tstart,tend,tlen,theader,gapopen,nident,mismatch",
        evalue="1.000E-06",
    log:
        outdir / "logs/homolog_prospection/candidates_filtering/mmseqs/{sample}.log",
    threads: 4
    conda:
        "../../envs/mmseqs2.yaml"
    shell:
        """
        mkdir -p temp_{wildcards.sample}
        mmseqs search {input.query} {input.probes} {wildcards.sample}_results temp_{wildcards.sample} --threads {threads} -e {params.evalue} --remove-tmp-files -a > {log} 2>&1
        mmseqs convertalis {input.query} {input.probes} {wildcards.sample}_results {output} --format-mode 4 --format-output {params.fields} --threads {threads} >> {log} 2>&1
        rm -r temp_{wildcards.sample}
        rm *_results.*
        """


rule candidates_filtering:
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
        min_cov=candidate_scf_min_cov,
        min_idt=candidate_scf_min_idt,
        separate_probes=False,
        tpat=pattern,
    script:
        "../../../src/homolog_filtering.py"
