from hmm_build import hmm_build

targets.append(
    expand(outdir / "translated_probes/probe_profiles/{probe}.hmm", probe=probes_list)
)


checkpoint split_probes:
    input:
        probes=outdir / "translated_probes/longest_cds.fasta",
    output:
        expand(
            outdir / "translated_probes/split_probes/{probe}.fasta", probe=probes_list
        ),
    run:
        split_dir = Path(output[0]).parent
        split_dir.mkdir(exist_ok=True, parents=True)
        for probe, outfile in zip(SeqIO.parse(input.probes, "fasta"), output):
            SeqIO.write(probe, outfile, "fasta")


rule build_probe_hmms:
    input:
        outdir / "translated_probes/split_probes/{probe}.fasta",
    output:
        outdir / "translated_probes/probe_profiles/{probe}.hmm",
    log:
        outdir / "logs/translated_probes/probe_profiles/{probe}.log",
    run:
        hmm_build(Path(input[0]), Path(output[0]), "amino")
