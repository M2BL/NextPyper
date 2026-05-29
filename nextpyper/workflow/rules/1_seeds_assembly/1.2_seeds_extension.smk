pathvars:
    stage="1_assembly/11_extension",
    results="<workdir>/<stage>",

rule make_diamond_probes_db:
    input:
        "<workdir>/0_preprocessing/longest_cds.fasta",
    output:
        "<workdir>/1_assembly/probes.dmnd",
    log:
        "<logs>/1_assembly/diamond_makedb.log",
    params:
        subpath(output[0], strip_suffix=".dmnd"),
    conda:
        "../../envs/matching.yaml"
    shell:
        "diamond makedb --db {params} --in {input} > {log} 2>&1"


rule raw_assembly_to_probes_matching:
    input:
        probes="<workdir>/1_assembly/probes.dmnd",
        query="<workdir>/1_assembly/10_raw_assembly/101_scaffolds/{sample}.fasta",
    output:
        "<results>/110_raw_matching_tables/{sample}.tsv",
    params:
        fields=lookup("diamond_matching/fields", within=pipeline),
        sensitivity=lookup("diamond_matching/sensitivity", within=pipeline),
        evalue=lookup("diamond_matching/evalue", within=pipeline),
        max_hsps=lookup("diamond_matching/max_hsps", within=pipeline),
        gapopen=lookup("diamond_matching/gapopen", within=pipeline),
        frameshift=lookup("diamond_matching/frameshift", within=pipeline),
        min_orf_len=lookup("diamond_matching/min_orf_len", within=pipeline),
    log:
        "<logs>/<stage>/110_raw_matching_diamond/{sample}.log",
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


def read_insert_size(stats: Path) -> int:
    insert_hist = json.loads(stats.read_text())["insert_size"]["histogram"]
    return last(filter(itemgetter(1), enumerate(insert_hist)))[0]


def get_max_intron_size(wildcards, input):
    """Parametrize the maximum intron size that can be bridged during extension
    by taking twice the maximum observed insert size of the data."""

    return read_insert_size(Path(input.stats)) * 2


rule extend_paths:
    input:
        graph="<workdir>/1_assembly/10_raw_assembly/100_spades/{sample}/assembly_graph_with_scaffolds.gfa",
        table="<results>/110_raw_matching_tables/{sample}.tsv",
        stats="<logs>/0_preprocessing/01_trimming_fastp/{sample}.json",
    output:
        ext_table="<results>/111_paths_extension/1110_tables/{sample}.tsv",
        ext_seqs="<results>/111_paths_extension/1111_sequences/{sample}.fasta",
    params:
        filter_low_dp_comps=lookup("filter_low_dp_comps", within=scf_ext),
        floor_len=lookup("floor_len_extension", within=scf_ext),
        ceil_len=lookup("ceiling_len_extension", within=scf_ext),
        plen_scaling=lookup("probe_len_scaling", within=scf_ext),
        max_extensions=lookup("max_extensions", within=scf_ext),
        min_idt=lookup("scf_min_idt", within=seeds_filt_params),
        max_intron_size=(
            get_max_intron_size
            if (size := lookup("max_intron_size", within=scf_ext)) == "auto"
            else size
        ),
        probe_pattern=lambda wildcards: probe_pattern,
    log:
        "<logs>/<stage>/111_path_extension/{sample}.log",
    script:
        "../../../src/gfa_graph.py"
