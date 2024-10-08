SAUTE_PATTERN = r"^Contig_(?P<sample>.*?)-(?P<probe>.*?)_(?P<cluster>\d+?)_(?P<seed>\d+?):(?P<component>\d+?):[^ ]+$"

BLAST_COLS = (
    "query",
    "subject",
    "idt",
    "aln_len",
    "mismatches",
    "gap_opens",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "score",
    "query_len",
    "matches",
)
BLAST_USECOLS = (
    "query",
    "subject",
    "mismatches",
    "gap_opens",
    "qstart",
    "qend",
    "query_len",
    "matches",
)


def orient_scf(rec: SeqRecord, trans: bool) -> SeqRecord:
    if trans:
        rec.seq = rec.seq.reverse_complement()
    return rec


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


rule make_blast_probe_db:
    input:
        outdir / "homolog_prospection/matching_probes.fasta",
    output:
        multiext(
            str(
                outdir / "homolog_prospection/blast_filtering/db/matching_probes.fasta"
            ),
            ".pdb",
            ".phr",
            ".pin",
            ".pot",
            ".psq",
            ".ptf",
            ".pto",
        ),
    params:
        outdir / "homolog_prospection/blast_filtering/db/matching_probes.fasta",
    log:
        outdir / "logs/homolog_prospection/blast_filtering/makedb.log",
    conda:
        "../../envs/blast.yaml"
    shell:
        "makeblastdb -dbtype prot -in {input} -out {params} 2> {log}"


rule scfs_blast_filtering:
    input:
        query=outdir / "saute/target_assembly/{sample}/target_vars.fasta",
        db=rules.make_blast_probe_db.output,
    output:
        outdir / "homolog_prospection/blast_filtering/blastx/{sample}.blast.tsv",
    params:
        db=outdir / "homolog_prospection/blast_filtering/db/matching_probes.fasta",
        others="-outfmt '6 std qlen nident'",
    log:
        outdir / "logs/homolog_prospection/blast_filtering/blastx/{sample}.log",
    threads: 4
    conda:
        "../../envs/blast.yaml"
    shell:
        "blastx -num_threads {threads} -query {input.query} -db {params.db} -out {output} {params.others} 2> {log}"


rule parse_blast_filtering:
    input:
        scfs=outdir / "saute/target_assembly/{sample}/target_vars.fasta",
        blast=outdir / "homolog_prospection/blast_filtering/blastx/{sample}.blast.tsv",
    output:
        outdir / "homolog_prospection/blast_filtering/filtered_scfs/{sample}.fasta",
    params:
        min_cov=homolog_scf_min_cov,
        min_idt=homolog_scf_min_idt,
    run:
        min_idt = params.min_idt
        min_cov = params.min_cov

        sel_cols = [
            "query",
            "subject",
            "cis",
            "matches",
            "mismatches",
            "gap_opens",
            "query_len",
        ]
        group_cols = ["query", "subject", "cis"]


        df = pl.read_csv(
            input.blast, separator="\t", has_header=False, new_columns=BLAST_COLS
        ).select(BLAST_USECOLS)

        final_scfs = (
            df.with_columns(
                qprobe=pl.col("query").str.extract(SAUTE_PATTERN, 2),
                cis=pl.col("qend") > pl.col("qstart"),
            )
            .filter(pl.col("subject").str.contains(pl.col("qprobe"), literal=True))
            .select(sel_cols)
            .group_by(group_cols)
            .agg(
                pl.sum("matches"),
                pl.sum("mismatches"),
                pl.sum("gap_opens"),
                pl.first("query_len"),
            )
            .with_columns(
                cov=(
                    (pl.col("matches") + pl.col("mismatches")) * 3 / pl.col("query_len")
                ),
                idt=(
                    pl.col("matches")
                    / (pl.col("matches") + pl.col("mismatches") + pl.col("gap_opens"))
                ),
            )
            .filter((pl.col("cov") > min_cov) & (pl.col("idt") > min_idt))
            .group_by("query")
            .agg(pl.all().sort_by("idt").last())
            .select(["query", "cis"])
            .with_columns(~pl.col("cis"))
        )

        filt_ids = dict(final_scfs.iter_rows())

        filt_scfs = (
            orient_scf(rec, trans)
            for rec in SeqIO.parse(input.scfs, "fasta")
            if (trans := filt_ids.get(rec.id)) is not None
        )
        SeqIO.write(filt_scfs, output[0], "fasta")
