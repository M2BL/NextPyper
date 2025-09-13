rule make_diamond_probes_db:
    input:
        outdir / "translated_probes/longest_cds.fasta",
    output:
        db=outdir / "assembled/filtering/dbs/probes.dmnd",
    log:
        outdir / "logs/assembled/filtering/dbs/diamond.log",
    params:
        db=subpath(output.db, strip_suffix=".dmnd"),
    conda:
        "../../envs/matching.yaml"
    shell:
        "diamond makedb --db {params.db} --in {input} > {log} 2>&1"


# rule make_mmseqs_probes_db:
#     input:
#         outdir / "translated_probes/longest_cds.fasta",
#     output:
#         outdir / "assembled/filtering/dbs/probes",
#     log:
#         outdir / "logs/assembled/filtering/dbs/probes.log",
#     conda:
#         "../../envs/matching.yaml"
#     shell:
#         "mmseqs createdb --dbtype 1 {input} {output} > {log} 2>&1 "


rule raw_assembly_to_probes_matching:
    input:
        probes=outdir / "assembled/filtering/dbs/probes.dmnd",
        # query=outdir / "assembled/extension/{sample}.fasta",
        query=outdir / "assembled/scaffolds/{sample}.fasta",
    output:
        # outdir / "assembled/filtering/matching_tables/{sample}.tsv",
        outdir / "assembled/filtering/raw_matching_tables/{sample}.tsv",
    params:
        fields=lookup("diamond_matching/fields", within=pipeline),
        sensitivity=lookup("diamond_matching/sensitivity", within=pipeline),
        evalue=lookup("diamond_matching/evalue", within=pipeline),
        max_hsps=lookup("diamond_matching/max_hsps", within=pipeline),
        gapopen=lookup("diamond_matching/gapopen", within=pipeline),
        frameshift=lookup("diamond_matching/frameshift", within=pipeline),
        min_orf_len=lookup("diamond_matching/min_orf_len", within=pipeline),
    log:
        # outdir / "logs/assembled/filtering/diamond/{sample}.log",
        outdir / "logs/assembled/filtering/raw_filtering/{sample}.log",
    threads: 4
    conda:
        "../../envs/matching.yaml"
    shell:
        """
        diamond blastx \
        --threads {threads} \
        --db {input.probes} --query {input.query} --out {output} \
        {params.sensitivity} \
        --evalue {params.evalue} \
        --max-hsps {params.max_hsps} \
        --min-orf {params.min_orf_len} \
        --gapopen {params.gapopen} \
        --frameshift {params.frameshift} \
        --outfmt 6 {params.fields} > {log} 2>&1
        """


# rule raw_assembly_to_probes_matching:
#     input:
#         probes=outdir / "assembled/filtering/dbs/probes",
#         query=outdir / "assembled/scaffolds/{sample}.fasta",
#     output:
#         outdir / "assembled/filtering/raw_matching_tables/{sample}.tsv",
#     params:
#         fields=lookup("mmseqs_matching/fields", within=pipeline),
#         evalue=lookup("mmseqs_matching/evalue", within=pipeline),
#         min_orf_len=lookup("mmseqs_matching/min_orf_len", within=pipeline),
#         sensitivity=lookup("mmseqs_matching/sensitivity", within=pipeline),
#     log:
#         outdir / "logs/assembled/filtering/raw_filtering/{sample}.log",
#     threads: 4
#     shadow:
#         "shallow"
#     conda:
#         "../../envs/matching.yaml"
#     shell:
#         """
#         mkdir -p temp_{wildcards.sample}
#         mmseqs easy-search {input.query} {input.probes} {output} temp_{wildcards.sample} \
#         --threads {threads} \
#         -s {params.sensitivity} \
#         -e {params.evalue} \
#         --min-length {params.min_orf_len} \
#         --format-mode 4 --format-output {params.fields} \
#         --remove-tmp-files -a > {log} 2>&1
#         """


def read_insert_size(stats: Path) -> int:
    insert_hist = json.loads(stats.read_text())["insert_size"]["histogram"]
    return last(filter(itemgetter(1), enumerate(insert_hist)))[0]


def get_max_intron_size(wildcards, input):
    """Parametrize the maximum intron size that can be bridged during extension
    by taking twice the maximum observed insert size of the data."""

    return read_insert_size(Path(input.stats)) * 2


rule extend_paths:
    input:
        graph=outdir / "assembled/spades/{sample}/assembly_graph_with_scaffolds.gfa",
        table=outdir / "assembled/filtering/raw_matching_tables/{sample}.tsv",
        stats=outdir / "logs/preprocessing/fastp/{sample}.json",
    output:
        outdir / "assembled/extension/{sample}.fasta",
    params:
        floor_len=lookup("floor_len_extension", within=scf_ext),
        ceil_len=lookup("ceiling_len_extension", within=scf_ext),
        plen_scaling=lookup("probe_len_scaling", within=scf_ext),
        max_extensions=lookup("max_extensions", within=scf_ext),
        max_intron_size=(
            get_max_intron_size
            if (size := lookup("max_intron_size", within=scf_ext)) == "auto"
            else size
        ),
        probe_pattern=lambda wildcards: probe_pattern,
    log:
        outdir / "logs/assembled/extension/{sample}.log",
    script:
        "../../../src/gfa_graph.py"
