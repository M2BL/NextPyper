rule prefix_scfs:
    input:
        branch(
            lookup(dpath="{sample}/type", within=sample_dict),
            cases={
                "rna": outdir / "assembled/spades/{sample}/transcripts.fasta",
                "targeted": outdir / "assembled/spades/{sample}/scaffolds.fasta",
                "": outdir / "assembled/spades/{sample}/scaffolds.fasta",
            },
        ),
    output:
        outdir / "assembled/prefixed/{sample}.fasta",
    run:
        prefix = f"{wildcards.sample}-"
        prefix_fasta(input[0], output[0], prefix)
