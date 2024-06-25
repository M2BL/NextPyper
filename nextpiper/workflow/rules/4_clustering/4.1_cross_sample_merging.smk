from pathlib import Path
from Bio import SeqIO
from contig_cluster import HDBcluster

graph_dir = outdir / Path(
    "assembled/simplyfied" if graph_simplification else "assembled/prefixed"
)


targets.append(expand(outdir / "clustering/clusters/{probes}", probes=probes_list))
# ToDo: Test:
# What happens if not all the probes are present in at least one sample?
# This would certainly happen in real life cases.


checkpoint split_probes:
    input:
        probes=outdir / "translated_probes/longest_cds.fasta",
    output:
        expand(outdir / "clustering/split_probes/{probe}.fasta", probe=probes_list),
    run:
        split_dir = Path(output[0]).parent
        split_dir.mkdir(exist_ok=True, parents=True)
        for probe, outfile in zip(SeqIO.parse(input.probes, "fasta"), output):
            SeqIO.write(probe, outfile, "fasta")


def union_probes(wildcards):
    glob_match = glob_wildcards(
        graph_dir / f"{{sample}}/{wildcards.probe}/contigs.fasta"
    )
    return expand(
        graph_dir / "{sample}/{probe}/contigs.fasta",
        sample=glob_match.sample,
        probe=wildcards.probe,
    )


rule merge_asms:
    input:
        direc=outdir / "logs/dones/prefixing.done",
        probes=union_probes,
    output:
        outfile=outdir / "clustering/sample_merged_input/{probe}.fasta",
    shell:
        "cat {input.probes} > {output}"


checkpoint clustering:
    input:
        probes=outdir / "clustering/split_probes/{probe}.fasta",
        contigs=outdir / "clustering/sample_merged_input/{probe}.fasta",
    output:
        directory(outdir / "clustering/clusters/{probe}"),
    log:
        outdir / "logs/clustering/{probe}.log",
    conda:
        "../../envs/clustering.yaml"
    script:
        "../../../src/contig_cluster.py"
