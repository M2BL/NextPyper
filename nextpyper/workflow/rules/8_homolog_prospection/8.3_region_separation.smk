POS_ALLELE_PATTERN = re.compile(
    r"^(?P<sample>.*?)-(?P<probe>.*?)_NODE_(?P<seed>\d+?):.*$", re.VERBOSE
)


def rename_rec(rec: SeqRecord, sample: str) -> SeqRecord:
    new_name = f"{sample}-{rec.id.split("-",1)[1]}"
    rec.id = rec.name = rec.description = new_name
    return rec


checkpoint per_probe_scaffold_grouping:
    input:
        expand(
            outdir / "homolog_prospection/allele_collapsing/vsearch/{sample}.fasta",
            sample=sample_list,
        ),
    output:
        directory(outdir / "homolog_prospection/region_separation/input_scfs"),
    log:
        outdir / "logs/homolog_prospection/region_separation/scfs_grouping.log",
    run:
        with open(log[0], "w") as outlog:
            sys.stdout = sys.stderr = outlog
            scfs = {Path(file).stem: list(SeqIO.parse(file, "fasta")) for file in input}
            renamed_scfs = {
                sample: [rename_rec(rec, sample) for rec in recs]
                for sample, recs in scfs.items()
            }
            all_recs = list(chain.from_iterable(renamed_scfs.values()))
            grouped_scfs = group_probes(
                all_recs, POS_ALLELE_PATTERN, match_group="probe"
            )

            outfolder = Path(output[0])
            outfolder.mkdir(exist_ok=True)
            for probe, recs in grouped_scfs.items():
                SeqIO.write(recs, outfolder / f"{probe}.fasta", "fasta")


checkpoint split_matching_probes:
    input:
        outdir / "homolog_prospection/matching_probes.fasta",
    output:
        directory(outdir / "homolog_prospection/region_separation/input_probes"),
    log:
        outdir / "logs/homolog_prospection/region_separation/probe_grouping.log",
    run:
        with open(log[0], "w") as outlog:
            sys.stdout = sys.stderr = outlog

            if multi_probes:
                probes = list(SeqIO.parse(input[0], "fasta"))
                outfolder = Path(output[0])
                outfolder.mkdir(exist_ok=True)
                for probe, recs in group_probes(probes, pattern).items():
                    SeqIO.write(recs, outfolder / f"{probe}.fasta", "fasta")

            else:
                outfolder = Path(output[0])
                outfolder.mkdir(exist_ok=True)
                ext_probe = re.compile(pattern)
                for probe_rec in SeqIO.parse(input[0], "fasta"):
                    probe = ext_probe.search(probe_rec.id)[1]
                    SeqIO.write(probe_rec, outfolder / f"{probe}.fasta", "fasta")


rule separate_cds_by_regions:
    input:
        probes_dir=outdir / "homolog_prospection/region_separation/input_probes",
        scfs_dir=outdir / "homolog_prospection/region_separation/input_scfs",
    output:
        directory(outdir / "homolog_prospection/region_separation/separation_output"),
    params:
        min_probe_contig_sim=min_probe_contig_sim,
        min_fragment_cov=min_fragment_cov,
        min_contig_length=min_contig_length,
    log:
        outdir / "logs/homolog_prospection/region_separation/separation.log",
    conda:
        "../../envs/clustering.yaml"
    script:
        "../../../src/miniprot.py"
