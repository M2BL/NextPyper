spades_pat = r"NODE_\d+_length_\d+_cov_(.*)"


rule prefilter_scfs_by_cov:
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
        outdir / "assembled/prefiltered/{sample}.fasta",
    params:
        min_cov=lookup("prefilter_by_cov/min_cov", within=pipeline),
    run:
        head_pat = re.compile(spades_pat)
        get_cov = lambda rec: float(head_pat.search(rec.id)[1])

        gen_recs = (
            rec
            for rec in SeqIO.parse(input[0], "fasta")
            if get_cov(rec) >= params.min_cov
        )
        SeqIO.write(gen_recs, output[0], "fasta")


rule prefix_scfs:
    input:
        outdir / "assembled/prefiltered/{sample}.fasta",
    output:
        outdir / "assembled/prefixed/{sample}.fasta",
    run:
        prefix = f"{wildcards.sample}-"
        prefix_fasta(input[0], output[0], prefix)
