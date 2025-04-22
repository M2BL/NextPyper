def select_asm_kmer(wildcards, input):
    if spades_k == "auto":
        with open(input.json) as file:
            summary = json.load(file)
            k_mid = int(summary["summary"]["after_filtering"]["read2_mean_length"]) // 2
            if k_mid >= 117:
                return "107,127"
            elif k_mid % 2 == 0:
                return f"{k_mid-9},{k_mid+11}"
            else:
                return f"{k_mid-10},{k_mid+10}"
    else:
        return spades_k


rule spades_assembly:
    input:
        in1=outdir / "preprocessed/trimmed/{sample}_R1.fastq.gz",
        in2=outdir / "preprocessed/trimmed/{sample}_R2.fastq.gz",
        json=outdir / "logs/preprocessing/fastp/{sample}.json",
    output:
        out_dir=directory(outdir / "assembled/spades/{sample}"),
        graph=outdir / "assembled/spades/{sample}/assembly_graph_with_scaffolds.gfa",
        contigs=outdir / "assembled/spades/{sample}/scaffolds.fasta",
    params:
        mode=lambda wildcards: (
            "--rna" if sample_dict[wildcards.sample]["type"] == "rna" else "--meta"
        ),
        k=select_asm_kmer,
        params=f"--only-assembler",
    log:
        outdir / "logs/assembled/spades/{sample}.log",
    threads: 4
    conda:
        "../../envs/assembly_spades.yaml"
    shell:
        "spades.py -t {threads} {params.mode} {params.params} -k {params.k} -1 {input.in1} -2 {input.in2} -o {output.out_dir} > {log} 2>&1"
