SAUTE_PATTERN = re.compile(
    r"^Contig_.*?-(?P<probe>.*?)_(?P<cluster>\d+?)_\d+?:\d+:[^ ]+$",
    re.VERBOSE,
)

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
        outdir / "log/homolog_prospection/blast_filtering/makedb.log",
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
        outdir / "log/homolog_prospection/blast_filtering/blastx/{sample}.log",
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
        filt_ids = {}

        df = pd.read_csv(input.blast, sep="\t", names=BLAST_COLS, usecols=BLAST_USECOLS)
        for query in df["query"].unique():
            probe = re.search(SAUTE_PATTERN, query)["probe"]
            dfq = df.query("query == @query")
            dfq = dfq[dfq.subject.str.contains(probe, regex=True)]

            if len(dfq) == 0:
                continue

            query_len = dfq["query_len"].to_list()[0]
            dfq["cis"] = dfq.eval("qend > qstart")
            agg = (
                dfq.loc[:, ["subject", "cis", "matches", "mismatches", "gap_opens"]]
                .groupby(by=["subject", "cis"])
                .sum()
            )
            metrics = pd.DataFrame(
                {
                    "cov": agg.eval("matches + mismatches") * 3 / query_len,
                    "idt": agg.eval("matches / (matches + mismatches + gap_opens)"),
                }
            )

            if len(result := metrics.query("idt >= @min_idt and cov >= @min_cov")) > 0:
                filt_ids[query] = not result.reset_index()["cis"][0]

        filt_scfs = (
            orient_scf(rec, trans)
            for rec in SeqIO.parse(input.scfs, "fasta")
            if (trans := filt_ids.get(rec.id)) is not None
        )
        SeqIO.write(filt_scfs, output[0], "fasta")
