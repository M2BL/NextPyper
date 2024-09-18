from gfa_parser import split_into_hmms, filter_components_hmm
from itertools import chain


def aggregate_hmms_spades(wildcards):
    if multi_probes:
        matches = {
            p: glob_wildcards(
                Path(checkpoints.make_probe_consensus.get(probe=p).output[0])
                / f"{p}_{{cluster}}.fasta"
            )
            for p in probes_list
        }

        return list(
            chain.from_iterable(
                expand(
                    outdir / f"translated_probes/probe_profiles/{p}_{{cluster}}.hmm",
                    cluster=matches[p].cluster,
                )
                for p in probes_list
            )
        )
    else:
        return expand(
            outdir / "translated_probes/probe_profiles/{probe}.hmm", probe=probes_list
        )


rule spades_assembly:
    input:
        in1=outdir / "preprocessed/filtered/{sample}_R1.fastq",
        in2=outdir / "preprocessed/filtered/{sample}_R2.fastq",
        hmms=aggregate_hmms_spades,
    output:
        out_dir=directory(outdir / "assembled/spades/{sample}"),
        contigs=outdir / "assembled/spades/{sample}/scaffolds.fasta",
        gfa=outdir / "assembled/spades/{sample}/assembly_graph_after_simplification.gfa",
        stats=outdir / "assembled/spades/{sample}/hmm_statistics.txt",
    params:
        params=f"--only-assembler --cov-cutoff auto {spades_k}",
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
    params:
        probes_dir=outdir / "translated_probes/multi_probe_consensus",
        min_probe_cov=0.1,
    output:
        directory(outdir / "assembled/split_components/{sample}"),
    run:
        probe_lens = {
            file.stem: len(SeqIO.read(file, "fasta")) * 3
            for file in Path(params.probes_dir).glob("*/*.fasta")
        }

        components = filter_components_hmm(
            Path(input.gfa), Path(input.hmm), probe_lens, params.min_probe_cov
        )
        split_into_hmms(
            gfa_path=Path(input.gfa),
            components=components,
            outdir=Path(output[0]),
            prefix=f"{wildcards.sample}-",
            write_graphs=False,
            write_seqs=True,
        )


## See Rule 3.1 for further explanation
def aggregate_split(wildcards):
    checkpoint_output = checkpoints.split_graph_into_hmms.get(
        sample=wildcards.sample
    ).output[0]
    global_match = glob_wildcards(Path(checkpoint_output) / "{probe}.fasta")
    return expand(
        outdir / f"assembled/split_components/{wildcards.sample}/{{probe}}.fasta",
        probe=global_match.probe,
    )


rule collect_prefixed_assemblies:
    input:
        aggregate_split,
    output:
        chkpt=outdir / "logs/assembled/collect_split/{sample}.chkpt",
    shell:
        "echo {input} | tr '[:space:]' '\n' >> {output.chkpt}"


checkpoint done_splitting:
    input:
        chkpt=expand(
            outdir / "logs/assembled/collect_split/{sample}.chkpt", sample=sample_list
        ),
    output:
        done=touch(outdir / "logs/dones/splitting.done"),
