from trans_dec_parser import transdecoder_parser

targets.append(outdir / "translated_probes/longest_cds.fasta")


rule transdecoder_longorfs:
    input:
        probes.resolve(),
    output:
        outdir / "translated_probes/transdecoder/longest_orfs.pep",
    params:
        out_folder=lambda wildcards, output: Path(output[0]).parent,
        min_len=100,  ## Todo: find a good compromise for this parameter
    log:
        outdir / "logs/translated_probes/transdecoder.log",
    conda:
        "../../envs/translating.yaml"
    shell:
        "TransDecoder.LongOrfs -m {params.min_len} -t {input} -O {params.out_folder} 2> {log}"


# ToDo: Review the current parsing, with too many orfs it is taking forever
rule parse_translation:
    input:
        pep=outdir / "translated_probes/transdecoder/longest_orfs.pep",
    output:
        outdir / "translated_probes/longest_cds.fasta",
    log:
        outdir / "logs/translated_probes/trans_dec_parser.log",
    threads: 1
    run:
        transdecoder_parser(str(input.pep), str(output))
