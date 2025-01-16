rule prefix_and_filter_scfs_by_cov:
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
    params:
        min_cov=min_scf_cov,
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """bioawk -c fastx '
        {{split($name, parts, "_"); 
        if((parts[6]*1)>=4)
        {{print ">{wildcards.sample}-"$name; $seq }} }}' {input} > {output} 
        """
