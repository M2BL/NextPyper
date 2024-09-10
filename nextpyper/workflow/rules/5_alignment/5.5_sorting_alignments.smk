from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def get_cluster(rec: SeqRecord) -> str:
    hit = saute_seq_pattern.match(rec.name)
    return f"{hit["probe"]}_{hit["cluster"]}"


def sort_aln_by_cluster(aln_path: Path, output_path: Path) -> None:
    aln = AlignIO.read(aln_path, "fasta")
    new_aln = MultipleSeqAlignment(sorted(aln, key=get_cluster))
    AlignIO.write(new_aln, output_path, "fasta")


rule sort_var_alns:
    input:
        outdir / "var_aligned/trimal/{probe}.fasta",
    output:
        outdir / "var_aligned/sorted_alns/{probe}.fasta",
    run:
        sort_aln_by_cluster(input[0], output[0])
