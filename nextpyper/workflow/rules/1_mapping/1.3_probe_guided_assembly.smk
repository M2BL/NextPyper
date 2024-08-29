from gfa_parser import split_into_hmms, filter_components_hmm

targets.append(
    expand(outdir / "assembled/split_components/{sample}", sample=sample_list)
)


rule spades_assembly:
    input:
        in1=outdir / "trimmed/{sample}_R1.fastq",
        in2=outdir / "trimmed/{sample}_R2.fastq",
        approv=outdir / "logs/dones/probe_hmms.done",
    output:
        out_dir=directory(outdir / "assembled/spades/{sample}"),
        contigs=outdir / "assembled/spades/{sample}/scaffolds.fasta",
        gfa=outdir / "assembled/spades/{sample}/assembly_graph_after_simplification.gfa",
        stats=outdir / "assembled/spades/{sample}/hmm_statistics.txt",
    params:
        params="--only-assembler --cov-cutoff auto",
        # params="--only-assembler --cov-cutoff auto -k 55,77", #ToDo: Consider for Speed-up
        hmms=outdir / "translated_probes/probe_profiles",
    log:
        outdir / "logs/assembled/spades/{sample}.log",
    threads: max(1, max_threads // len(sample_list))
    conda:
        "../../envs/assembly_spades.yaml"
    shell:
        "spades.py -t {threads} {params.params} -1 {input.in1} -2 {input.in2} --custom-hmms {params.hmms} -o {output.out_dir} > {log} 2>&1"


checkpoint split_graph_into_hmms:
    input:
        gfa=outdir / "assembled/spades/{sample}/assembly_graph_after_simplification.gfa",
        hmm=outdir / "assembled/spades/{sample}/hmm_statistics.txt",
    output:
        directory(outdir / "assembled/split_components/{sample}"),
    run:
        components = filter_components_hmm(Path(input.gfa), Path(input.hmm))
        split_into_hmms(Path(input.gfa), components, Path(output[0]))
