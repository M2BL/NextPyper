def select_asm_kmer(wildcards, input):
    spades_k = lookup("spades/k", within=pipeline)
    if spades_k == "auto":
        summary = json.loads(Path(input.json).read_text())
        k_mid = int(summary["summary"]["after_filtering"]["read2_mean_length"]) // 2
        if k_mid >= 117:
            return "21,45,107,127"
        elif k_mid % 2 == 0:
            return f"21,45,{k_mid-9},{k_mid+11}"
        else:
            return f"21,45,{k_mid-10},{k_mid+10}"
    else:
        return spades_k


rule spades_assembly:
    input:
        in1=outdir / "preprocessed/cleaned/{sample}_R1.fastq.gz",
        in2=outdir / "preprocessed/cleaned/{sample}_R2.fastq.gz",
        json=outdir / "logs/preprocessing/fastp/{sample}.json",
    output:
        graph=outdir / "assembled/spades/{sample}/assembly_graph_with_scaffolds.gfa",
    params:
        k=select_asm_kmer,
        out_dir=subpath(output.graph, parent=True),
        extra="--meta --only-assembler",
    log:
        outdir / "logs/assembled/spades/{sample}.log",
    threads: 4
    shadow:
        "minimal"
    conda:
        "../../envs/assembly_spades.yaml"
    shell:
        "spades.py -t {threads} {params.extra} -k {params.k} -1 {input.in1} -2 {input.in2} -o {params.out_dir} > {log} 2>&1"


rule make_assembly_scaffolds:
    input:
        outdir / "assembled/spades/{sample}/assembly_graph_with_scaffolds.gfa",
    output:
        outdir / "assembled/scaffolds/{sample}.fasta",
    script:
        "../../../src/gfa2fasta.py"
