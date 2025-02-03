rule make_mmseqs_raw_assembly_dbs:
    input:
        branch(
            lookup(dpath="{sample}/type", within=sample_dict),
            cases={
                "rna": outdir / "assembled/spades/{sample}/transcripts.fasta",
                "targeted": outdir / "assembled/spades/{sample}/scaffolds.fasta",
                "": outdir / "assembled/spades/{sample}/scaffolds.fasta",
            },
        ),
    output:
        outdir / "assembled/filtering/dbs/raw_assembly/{sample}",
    log:
        outdir / "logs/assembled/filtering/make_raw_assembly_db/{sample}.log",
    conda:
        "../../envs/mmseqs2.yaml"
    shell:
        "mmseqs createdb --dbtype 2 {input} {output} > {log} 2>&1"


rule raw_assembly_to_probes_matching:
    input:
        probes=outdir / "assembled/filtering/dbs/probes/probes",
        query=outdir / "assembled/filtering/dbs/raw_assembly/{sample}",
    output:
        outdir / "assembled/filtering/raw_matching_tables/{sample}.tsv",
    params:
        fields=mmseq_fields,
        evalue=mmseq_evalue,
        min_orf_len=min_orf_len,
        sensitivity=mmseq_prefilt_sens,
    log:
        outdir / "logs/assembled/filtering/raw_filtering/{sample}.log",
    threads: 4
    conda:
        "../../envs/mmseqs2.yaml"
    shell:
        """
        mkdir -p temp_{wildcards.sample}
        mmseqs search {input.query} {input.probes} {wildcards.sample}_results temp_{wildcards.sample} --threads {threads} -s {params.sensitivity} -e {params.evalue} --min-length {params.min_orf_len} --remove-tmp-files -a > {log} 2>&1
        mmseqs convertalis {input.query} {input.probes} {wildcards.sample}_results {output} --format-mode 4 --format-output {params.fields} --threads {threads} >> {log} 2>&1
        rm -r temp_{wildcards.sample}
        rm {wildcards.sample}_results.*
        """


rule extend_paths:
    input:
        graph=outdir / "assembled/spades/{sample}/assembly_graph_with_scaffolds.gfa",
        table=outdir / "assembled/filtering/raw_matching_tables/{sample}.tsv",
    output:
        outdir / "assembled/extension/{sample}.fasta",
    params:
        floor_len=floor_len_extension,
        plen_scaling=plen_scaling_factor,
    log:
        outdir / "logs/assembled/extension/{sample}.log",
    script:
        "../../../src/gfa_graph.py"


rule prefix_and_filter_scfs_by_cov:
    input:
        outdir / "assembled/extension/{sample}.fasta",
    output:
        outdir / "assembled/prefixed/{sample}.fasta",
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """bioawk -c fastx '{{split($name, parts, "_"); 
        printf ">{wildcards.sample}-"; 
        for(i=1; i<=length(parts)-2; i++) 
        {{printf "%s%s", (i>1?"_":""), parts[i]}}; 
        print "\\n"$seq }}' {input} > {output}
        """
