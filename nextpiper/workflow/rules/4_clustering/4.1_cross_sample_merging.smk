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


def union_probes(wildcards):
    probes = {file.parent.name for file in graph_dir.glob("*/*/contigs.fasta")}
    probe_dict = {
        probe: list(graph_dir.glob(f"*/{probe}/contigs.fasta")) for probe in probes
    }
    return probe_dict[wildcards.probe]


checkpoint merge_asms:
    input:
        unpack(union_probes),
    output:
        outdir / "clustering/sample_merged_input/{probe}.fasta",
    shell:
        "cat {input} > {output}"


checkpoint split_probes:
    input:
        probes=outdir / "translated_probes/longest_cds.fasta",
    output:
        directory(outdir / "clustering/split_probes"),
    run:
        split_dir = Path(output[0])
        split_dir.mkdir(exist_ok=True, parents=True)
        for probe in SeqIO.parse(input.probes, "fasta"):
            SeqIO.write(probe, split_dir / f"{probe.name}.fasta", "fasta")


def aggregate_probes(wildcards):
    checkpoint_output = checkpoints.split_probes.get(**wildcards).output[0]
    return expand(
        outdir / "clustering/split_probes/{probe}.fasta",
        probe=glob_wildcards(os.path.join(checkpoint_output, "{probe}.fasta")).probe,
    )


checkpoint clustering:
    input:
        probes=outdir / "clustering/split_probes/{probe}.fasta",
        # probes=aggregate_probes,
        contigs=outdir / "clustering/sample_merged_input/{probe}.fasta",
    output:
        directory(outdir / "clustering/clusters/{probe}"),
    log:
        outdir / "logs/clustering/{probe}.log",
    conda:
        "../../envs/clustering.yaml"
    script:
        "../../../src/contig_cluster.py"
