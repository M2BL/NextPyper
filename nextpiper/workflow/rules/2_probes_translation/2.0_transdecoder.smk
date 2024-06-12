from trans_dec_parser import transdecoder_parser

targets.append(outdir / "translated_probes/longest_cds.fasta")


rule transdecoder_longorfs:
    input:
        outdir / "mapped/map_index/probes.fasta",
    output:
        directory(outdir / "translated_probes/transdecoder"),
    log:
        outdir / "logs/translated_probes/transdecoder.log",
    conda:
        "../../envs/translating.yaml"
    shell:
        "TransDecoder.LongOrfs -t {input} -O {output} 2> {log}"


rule parse_translation:
    input:
        pep=outdir / "translated_probes/transdecoder",
    output:
        outdir / "translated_probes/longest_cds.fasta",
    log:
        outdir / "logs/translated_probes/trans_dec_parser.log",
    threads: 1
    run:
        input_file = str(input.pep) + "/longest_orfs.pep"
        transdecoder_parser(input_file, str(output))
