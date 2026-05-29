pathvars:
    stage="1_assembly/10_raw_assembly",
    results="<workdir>/<stage>",

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
        in1="<workdir>/0_preprocessing/02_cleaning/{sample}_R1.fastq.gz",
        in2="<workdir>/0_preprocessing/02_cleaning/{sample}_R2.fastq.gz",
        json="<logs>/0_preprocessing/01_trimming_fastp/{sample}.json",
    output:
        graph="<results>/100_spades/{sample}/assembly_graph_with_scaffolds.gfa",
    params:
        k=select_asm_kmer,
        out_dir=subpath(output.graph, parent=True),
        extra="--meta --only-assembler",
    log:
        "<logs>/<stage>/100_spades/{sample}.log",
    threads: 4
    shadow:
        "minimal"
    conda:
        "../../envs/assembly.yaml"
    shell:
        "spades.py -t {threads} {params.extra} -k {params.k} -1 {input.in1} -2 {input.in2} -o {params.out_dir} > {log} 2>&1"


rule make_assembly_scaffolds:
    input:
        "<results>/100_spades/{sample}/assembly_graph_with_scaffolds.gfa",
    output:
        "<results>/101_scaffolds/{sample}.fasta",
    script:
        "../../../src/gfa2fasta.py"
