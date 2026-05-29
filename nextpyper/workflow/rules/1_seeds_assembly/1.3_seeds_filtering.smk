pathvars:
    stage="1_assembly/12_filtering",
    results="<workdir>/<stage>",
    inter="1_assembly/11_extension/111_paths_extension/1111_sequences"

use rule raw_assembly_to_probes_matching as seeds_to_probes_matching with:
    input:
        probes="<workdir>/1_assembly/probes.dmnd",
        query="<workdir>/<inter>/{sample}.fasta",
    output:
        "<results>/120_matching_tables/{sample}.tsv",
    log:
        "<logs>/<stage>/120_matching_diamond/{sample}.log",


rule seeds_coverage:
    input:
        scfs="<workdir>/<inter>/{sample}.fasta",
        clean1="<workdir>/0_preprocessing/02_cleaning/{sample}_R1.fastq.gz",
        clean2="<workdir>/0_preprocessing/02_cleaning/{sample}_R2.fastq.gz",
    output:
        counts="<results>/121_coverage/{sample}.counts",
        metabat="<results>/121_coverage/{sample}.metabat",
    params:
        extra="--proper-pairs-only --exclude-supplementary",
    log:
        "<logs>/<stage>/121_coverage/{sample}.log",
    threads: 4
    shadow:
        "shallow"
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """
        minimap2 -t {threads} -ax sr {input.scfs} {input.clean1} {input.clean2} 2> {log} | \
        samtools sort -u -@ {threads} > tmp_{wildcards.sample}.bam 2>> {log}
        coverm contig -m count {params.extra} -b tmp_{wildcards.sample}.bam  > {output.counts} 2>> {log}
        coverm contig -m metabat {params.extra} -b tmp_{wildcards.sample}.bam  > {output.metabat} 2>> {log}
        """


rule seeds_filtering:
    input:
        scfs="<workdir>/<inter>/{sample}.fasta",
        hits="<results>/120_matching_tables/{sample}.tsv",
        covs="<results>/121_coverage/{sample}.metabat",
    output:
        metrics="<results>/122_filt_seeds/1220_summary_tables/{sample}.tsv",
        scfs="<results>/122_filt_seeds/1221_filtered_scfs/{sample}.fasta",
    log:
        "<logs>/<stage>/122_seeds_filtering/{sample}.log",
    params:
        min_cov=lookup("scf_min_cov", within=seeds_filt_params),
        min_idt=lookup("scf_min_idt", within=seeds_filt_params),
        max_gc=lookup("max_gc", within=seeds_filt_params),
        min_gc=lookup("min_gc", within=seeds_filt_params),
        cov_threshold=lookup("cov_threshold", within=seeds_filt_params),
        cov_dynamic_filt=lookup("cov_dynamic_filt", within=seeds_filt_params),
        separate_probes=lambda wildcards: False,
        tag_scfs=lambda wildcards: False,
        qpat=lambda wildcards: SEED_PAT,
        tpat=lambda wildcards: probe_pattern,
    script:
        "../../../src/homolog_filtering.py"
