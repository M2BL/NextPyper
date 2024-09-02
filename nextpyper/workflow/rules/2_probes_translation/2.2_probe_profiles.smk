from hmm_build import hmm_build, hmm_consensus

targets.append(outdir / "logs/dones/probe_hmms.done")


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
            grouped_probes = group_probes(recs, pattern)
            out_dict = {Path(out).stem: out for out in output}

            for probe, recs in probe_recs.items():
                SeqIO.write(recs, out_dict[probe], "fasta")


rule mafft_probes:
    input:
        outdir / "translated_probes/split_probes/{probe}.fasta",
    output:
        outdir / "translated_probes/aln_probes/{probe}.fasta",
    params:
        "--auto ",
    log:
        outdir / "logs/translates_probes/aln_probes/{probe}.log",
    conda:
        "../../envs/alignment.yaml"
    shell:
        "mafft {params} {input} > {output} 2> {log}"


rule build_probe_hmms:
    input:
        branch(
            multi_probes,
            then=outdir / "translated_probes/aln_probes/{probe}.fasta",
            otherwise=outdir / "translated_probes/split_probes/{probe}.fasta",
        ),
    output:
        outdir / "translated_probes/probe_profiles/{probe}.hmm",
    log:
        outdir / "logs/translated_probes/probe_profiles/{probe}.log",
    run:
        hmm_build(Path(input[0]), Path(output[0]), "amino")


if multi_probes:

    rule build_probe_consensus:
        input:
            outdir / "translated_probes/probe_profiles/{probe}.hmm",
        output:
            consensus=outdir / "translated_probes/probe_consensus/{probe}.fasta",
        log:
            outdir / "logs/translated_probes/probe_consensus/{probe}.log",
        run:
            hmm_consensus(Path(input[0]), Path(output[0]))


rule done_probe_hmms:
    input:
        chkpt=expand(
            outdir / "translated_probes/probe_profiles/{probe}.hmm", probe=probes_list
        ),
    output:
        done=touch(outdir / "logs/dones/probe_hmms.done"),
