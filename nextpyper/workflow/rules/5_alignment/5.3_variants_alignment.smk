from collections import defaultdict
import re
from prefix_seqs import pref_rec

targets.append(outdir / "logs/dones/var_alns.done")


checkpoint group_variants_by_probe:
    input:
        expand(
            outdir / "saute/target_assembly/{samples}/target_vars.fasta",
            samples=sample_list,
        ),
    output:
        directory(outdir / "var_aligned/aln_input"),
    run:
        pattern = re.compile(r"Contig_(?P<probe>\w+)_(?P<cluster>\d+)")
        probe_recs = defaultdict(list)

        for sample_path in input:
            sample_path = Path(sample_path)
            for rec in SeqIO.parse(sample_path, "fasta"):
                probe = pattern.match(rec.name)["probe"]
                rec.name = rec.name.replace(":", "-")
                probe_recs[probe].append(pref_rec(rec, f"{sample_path.parent.stem}_"))

        print(probe_recs)
        output_dir = Path(output[0])
        output_dir.mkdir(parents=True, exist_ok=True)
        for probe, recs in probe_recs.items():
            SeqIO.write(recs, output_dir / f"{probe}.fasta", "fasta")


def get_vars_for_aln(wildcards):
    checkpoint_output = checkpoints.group_variants_by_probe.get(**wildcards).output[0]
    return Path(checkpoint_output) / f"{wildcards.probe}.fasta"


use rule mafft as mafft_vars with:
    input:
        alns=get_vars_for_aln,
    output:
        alns=outdir / "var_aligned/var_alns/{probe}.fasta",
    params:
        "--auto ",
    log:
        outdir / "logs/var_aligned/var_alns/{probe}.log",


# def aggregate_var_alns(wildcards):
#     checkpoint_output = checkpoints.group_variants_by_probe.get(
#         probe=wildcards.probe
#     ).output[0]
#     glob_match = glob_wildcards(Path(checkpoint_output) / f"{wildcards.probe}.fasta")

#     return expand(
#         outdir / "var_aligned/var_alns/{probe}.fasta",
#         probe=wildcards.probe,
#     )


## Is this collections necessary?
rule collect_var_alignments:
    input:
        # aggregate_var_alns,
        outdir / "var_aligned/var_alns/{probe}.fasta",
    output:
        chkpt=outdir / "logs/var_aligned/collect_var_alns/{probe}.chkpt",
    shell:
        "echo {input} | tr '[:space:]' '\n' >> {output.chkpt}"


checkpoint done_var_alns:
    input:
        chkpt=expand(
            outdir / "logs/var_aligned/collect_var_alns/{probe}.chkpt",
            probe=probes_list,
        ),
    output:
        done=touch(outdir / "logs/dones/var_alns.done"),
