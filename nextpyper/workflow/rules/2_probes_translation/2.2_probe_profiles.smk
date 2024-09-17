from hmm_build import hmm_build, hmm_consensus
from mmseqs import generate_consensuses


rule split_probes:
    input:
        probes=outdir / "translated_probes/longest_cds.fasta",
    output:
        expand(
            outdir / "translated_probes/split_probes/{probe}.fasta", probe=probes_list
        ),
    run:
        split_dir = Path(output[0]).parent
        split_dir.mkdir(exist_ok=True, parents=True)

        if not multi_probes:
            for probe, outfile in zip(SeqIO.parse(input.probes, "fasta"), output):
                SeqIO.write(probe, outfile, "fasta")
        else:
            probe_recs = list(SeqIO.parse(input.probes, "fasta"))
            grouped_probes = group_probes(probe_recs, pattern)
            out_dict = {Path(out).stem: out for out in output}

            for probe, recs in grouped_probes.items():
                SeqIO.write(recs, out_dict[probe], "fasta")


if multi_probes:

    rule mmseqs2_probe_consensus:
        input:
            probes=outdir / "translated_probes/split_probes/{probe}.fasta",
        output:
            outdir / "translated_probes/grouped_probes/{probe}_all_seqs.fasta",
        params:
            prefix=lambda wildcards: outdir
            / f"translated_probes/grouped_probes/{wildcards.probe}",
            other=f"--min-seq-id {mmseq2_min_seq_id}",
        threads: 1
        conda:
            "../../envs/mmseqs2.yaml"
        log:
            outdir / "logs/translated_probes/mmseq2/{probe}.log",
        shell:
            """
            mkdir -p tmp_mmseq_{wildcards.probe}
            mmseqs easy-linclust --threads {threads} {input.probes} {params.prefix} tmp_mmseq_{wildcards.probe} {params.other} > {log}
            rm -r tmp_mmseq_{wildcards.probe}
            """

    checkpoint make_probe_consensus:
        input:
            outdir / "translated_probes/grouped_probes/{probe}_all_seqs.fasta",
        output:
            directory(outdir / "translated_probes/multi_probe_consensus/{probe}"),
        run:
            file = Path(input[0])
            folder = Path(output[0])
            folder.mkdir(exist_ok=True, parents=True)
            for i, rec in enumerate(generate_consensuses(file)):
                name = f"{file.stem.removesuffix("_all_seqs")}_{i}.fasta"
                SeqIO.write(rec, folder / name, "fasta")

    def get_multi_probe_consensus(wildcards):
        w = wildcards
        chk_output = checkpoints.make_probe_consensus.get(probe=w.probe).output[0]
        match = glob_wildcards(Path(chk_output) / f"{w.probe}_{{cluster}}.fasta")

        return expand(
            outdir
            / f"translated_probes/multi_probe_consensus/{w.probe}/{w.probe}_{{cluster}}.fasta",
            cluster=match.cluster,
        )

    rule build_probe_hmms_multi:
        input:
            get_multi_probe_consensus,
        output:
            outdir / "translated_probes/probe_profiles/{probe}_{cluster}.hmm",
        # log:
        #     outdir / "logs/translated_probes/probe_profiles/{probe}_{cluster}.log",
        run:
            hmm_build(Path(input[0]), Path(output[0]), "amino")

    rule build_probe_consensus:
        input:
            outdir / "translated_probes/probe_profiles/{probe}_{cluster}.hmm",
        output:
            consensus=outdir
            / "translated_probes/probe_consensus/{probe}_{cluster}.fasta",
        # log:
        #     outdir / "logs/translated_probes/probe_consensus/{probe}_{cluster}.log",
        run:
            hmm_consensus(Path(input[0]), Path(output[0]))

else:

    rule build_probe_hmms_single:
        input:
            outdir / "translated_probes/split_probes/{probe}.fasta",
        output:
            outdir / "translated_probes/probe_profiles/{probe}.hmm",
        # log:
        #     outdir / "logs/translated_probes/probe_profiles/{probe}.log",
        run:
            hmm_build(Path(input[0]), Path(output[0]), "amino")


def aggregate_hmms(wildcards):
    w = wildcards
    if not multi_probes:
        return expand(
            outdir / f"translated_probes/probe_profiles/{w.probe}.hmm",
            probe=probes_list,
        )
    else:
        chk_output = checkpoints.make_probe_consensus.get(probe=w.probe).output[0]
        match = glob_wildcards(
            outdir / f"translated_probes/probe_profiles/{w.probe}_{{cluster}}.hmm"
        )

        return expand(
            outdir / f"translated_probes/probe_profiles/{w.probe}_{{cluster}}.hmm",
            cluster=match.cluster,
        )


rule collect_hmms:
    input:
        aggregate_hmms,
    output:
        chkpt=outdir / "logs/translated_probes/probe_profiles/{probe}.chkpt",
    shell:
        "touch {output.chkpt}"


rule done_probe_hmms:
    input:
        expand(
            outdir / "logs/translated_probes/probe_profiles/{probe}.chkpt",
            probe=probes_list,
        ),
    output:
        done=touch(outdir / "logs/dones/probe_hmms.done"),
