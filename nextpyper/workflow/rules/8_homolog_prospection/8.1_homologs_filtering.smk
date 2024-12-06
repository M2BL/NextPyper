rule gather_matching_probes:
    input:
        direc=outdir / "logs/dones/splitting.done",
        probes=outdir / "translated_probes/longest_cds.fasta",
    output:
        outdir / "homolog_prospection/matching_probes.fasta",
    params:
        split_dir=outdir / "assembled/split_components",
        probe_clusters_dir=outdir / "translated_probes/grouped_probes",
    run:
        subset_probes = set()
        cols = ["cluster", "probe"]

        probe_clusters = {
            file.stem for file in Path(params.split_dir).glob("*/*.fasta")
        }
        for probe_cluster in probe_clusters:
            probe, cluster = probe_cluster.rsplit("_", 1)
            clusters_path = f"{params.probe_clusters_dir}/{probe}_cluster.tsv"
            df = pd.read_csv(clusters_path, sep="\t", names=cols)
            probe_cluster = df.cluster.unique()[int(cluster)]
            subset_probes.update(set(df.query("cluster == @probe_cluster")["probe"]))

        probes_gen = (
            rec for rec in SeqIO.parse(input.probes, "fasta") if rec.id in subset_probes
        )
        SeqIO.write(probes_gen, Path(output[0]), "fasta")


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
        "mmseqs createdb --dbtype 1 {input} {output} 2> {log}"


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
        "mmseqs createdb --dbtype 2 {input} {output} 2> {log}"


rule candidates_to_probes_matching:
    input:
        probes=outdir
        / "homolog_prospection/candidates_filtering/dbs/probes/matching_probes",
        query=outdir / "homolog_prospection/candidates_filtering/dbs/samples/{sample}",
    output:
        outdir / "homolog_prospection/candidates_filtering/matching_tables/{sample}.tsv",
    params:
        db=outdir / "homolog_prospection/blast_filtering/db/matching_probes.fasta",
        others="-outfmt '6 std qlen nident'",
    log:
        outdir / "logs/homolog_prospection/blast_filtering/blastx/{sample}.log",
    threads: 4
    conda:
        "../../envs/mmseqs2.yaml"
    shell:
        """
        mkdir -p temp_{wildcards.sample}
        mmseqs search {input.query} {input.probes} {wildcards.sample}_results temp_{wildcards.sample} --threads {threads} -e 1.000E-06 --remove-tmp-files -a
        mmseqs convertalis {input.query} {input.probes} {wildcards.sample}_results {output} --format-mode 4 --format-output query,evalue,qstart,qend,qlen,tstart,tend,tlen,theader,gapopen,nident,mismatch --threads {threads}
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
    params:
        min_cov=homolog_scf_min_cov,
        min_idt=homolog_scf_min_idt,
    script:
        "../../../src/homolog_filtering.py"
