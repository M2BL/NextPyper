rule make_mmseqs_probes_db:
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


rule raw_assembly_to_probes_matching:
    input:
        probes=outdir / "assembled/filtering/dbs/probes/probes",
        query=outdir / "assembled/scaffolds/{sample}.fasta",
    output:
        outdir / "assembled/filtering/raw_matching_tables/{sample}.tsv",
    params:
        fields=mmseq_fields,
        evalue=mmseq_evalue,
        min_orf_len=min_orf_len,
        sensitivity=mmseq_sens,
    log:
        outdir / "logs/assembled/filtering/raw_filtering/{sample}.log",
    threads: 4
    conda:
        "../../envs/mmseqs2.yaml"
    shell:
        """
        mkdir -p temp_{wildcards.sample}
        mmseqs easy-search {input.query} {input.probes} {output} temp_{wildcards.sample} --threads {threads} -s {params.sensitivity} -e {params.evalue} --min-length {params.min_orf_len} --format-mode 4 --format-output {params.fields} --remove-tmp-files -a > {log} 2>&1
        rm -r temp_{wildcards.sample}
        """


def get_max_intron_size(wildcards, input):
    """Parametrize the maximum intron size that can be bridged during extension
    by taking twice the maximum observed insert size of the data."""

    with open(input.stats) as file:
        insert_hist = json.load(file)["insert_size"]["histogram"]
        return last(filter(itemgetter(1), enumerate(insert_hist)))[0] * 2


rule extend_paths:
    input:
        graph=outdir / "assembled/spades/{sample}/assembly_graph_with_scaffolds.gfa",
        table=outdir / "assembled/filtering/raw_matching_tables/{sample}.tsv",
        stats=outdir / "logs/preprocessing/fastp/{sample}.json",
    output:
        outdir / "assembled/extension/{sample}.fasta",
    params:
        floor_len=floor_len_extension,
        ceil_len=ceil_len_extension,
        plen_scaling=plen_scaling_factor,
        max_extensions=max_extensions,
        max_intron_size=(
            max_intron_size if max_intron_size != "auto" else get_max_intron_size
        ),
        probe_pattern=lambda wildcards: pattern,
    log:
        outdir / "logs/assembled/extension/{sample}.log",
    script:
        "../../../src/gfa_graph.py"


rule prefix_seeds:
    input:
        outdir / "assembled/extension/{sample}.fasta",
    output:
        outdir / "assembled/prefixed/{sample}.fasta",
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """bioawk -c fastx '{{printf ">{wildcards.sample}-%s",$name; print "\\n"$seq}}' {input} > {output}"""
