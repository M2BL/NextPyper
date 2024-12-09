rule spades_assembly:
    input:
        in1=outdir / "preprocessed/cleaned/{sample}_R1.fastq.gz",
        in2=outdir / "preprocessed/cleaned/{sample}_R2.fastq.gz",
    output:
        out_dir=directory(outdir / "assembled/spades/{sample}/"),
        contigs=outdir / "assembled/spades/{sample}/scaffolds.fasta",
    params:
        mode=lambda wildcards: (
            "--rna" if sample_dict[wildcards.sample]["type"] == "rna" else "--meta"
        ),
        params=f"--only-assembler -k {spades_k}",
    log:
        outdir / "logs/assembled/spades/{sample}.log",
    threads: max(1, max_threads // len(sample_list))
    conda:
        "../../envs/assembly_spades.yaml"
    shell:
        "spades.py -t {threads} {params.mode} {params.params} -1 {input.in1} -2 {input.in2} -o {output.out_dir} > {log} 2>&1"
