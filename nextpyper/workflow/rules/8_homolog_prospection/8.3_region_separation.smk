POS_ALLELE_PATTERN = re.compile(
    r"^(?P<sample>.*?)\|(?P<seed>.*?)-(?P<probe>.*?)_EDGE_(?P<seed_id>\d+)_length_(?P<len>\d+):[^ ]+$",
    re.VERBOSE,
)


def rename_rec(rec: SeqRecord, sample: str) -> SeqRecord:
    new_name = f"{sample}|{rec.id.removeprefix("Contig_")}"
    rec.id = rec.name = rec.description = new_name
    return rec


checkpoint per_probe_scaffold_grouping:
    input:
        expand(
            outdir / "homolog_prospection/allele_collapsing/vsearch/{sample}.fasta",
            sample=sample_list,
        ),
    output:
        expand(
            outdir / "homolog_prospection/region_separation/input_scfs/{probe}.fasta",
            probe=probes_list,
        ),
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

            outfolder = Path(output[0]).parent
            outfolder.mkdir(exist_ok=True)
            for probe, recs in grouped_scfs.items():
                SeqIO.write(recs, outfolder / f"{probe}.fasta", "fasta")

            for probe in probes_list:
                (outfolder / f"{probe}.fasta").touch(exist_ok=True)


checkpoint split_matching_probes:
    input:
        probes=outdir / "homolog_prospection/matching_probes.fasta",
        tables=expand(
            outdir
            / "logs/homolog_prospection/candidates_filtering/scfs_filtering/{samples}.log",
            samples=sample_list,
        ),
    output:
        expand(
            outdir
            / "homolog_prospection/region_separation/input_probes/{probe}.fasta",
            probe=probes_list,
        ),
    log:
        outdir / "logs/homolog_prospection/region_separation/probe_grouping.log",
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """
        cat {input.tables} | cut -f 2 | sort | uniq > probe_ids.txt
        seqkit grep -nf probe_ids.txt {input.probes} > temp_matching_probes.fasta 2> {log}
        rm probe_ids.txt

        for outfile in {output}; do
            name=$(basename $outfile .fasta)
            seqkit grep -rnp "$name" temp_matching_probes.fasta > $outfile
        done
        rm temp_matching_probes.fasta
        """


checkpoint separate_cds_by_regions:
    input:
        probes=outdir
        / "homolog_prospection/region_separation/input_probes/{probe}.fasta",
        scfs=outdir / "homolog_prospection/region_separation/input_scfs/{probe}.fasta",
    output:
        directory(
            outdir
            / "homolog_prospection/region_separation/separation_output/scfs/{probe}"
        ),
    params:
        min_probe_scaffold_sim=min_probe_scaffold_sim,
        min_fragment_cov=min_fragment_cov,
        min_exonic_length=min_exonic_length,
        substitution_matrix=blosum62,
    log:
        outdir / "logs/homolog_prospection/region_separation/separation/{probe}.log",
    conda:
        "../../envs/clustering.yaml"
    script:
        "../../../src/miniprot.py"


rule align_regions:
    input:
        outdir / "homolog_prospection/region_separation/separation_output/scfs/{probe}",
    output:
        directory(
            outdir / "homolog_prospection/region_separation/alns/{probe}",
        ),
    params:
        "--auto --adjustdirection",
    log:
        outdir / "logs/homolog_prospection/region_separation/alns/{probe}.log",
    threads: 1
    conda:
        "../../envs/alignment.yaml"
    shell:
        """
        rm -f {log}
        mkdir -p {output}

        for file in $(find {input} -name "*.fasta"); do
            name=$(basename $file)
            mafft {params} $file > {output}/$name 2>> {log}
        done
        """
