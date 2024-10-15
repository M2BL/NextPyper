from trans_dec_parser import transdecoder_parser


rule orffinder:
    input:
        probes_path.resolve(),
    output:
        outdir / "translated_probes/orfipy/longest_orfs.pep",
    params:
        out_folder=lambda wildcards, output: Path(output[0]).parent,
        pep=lambda wildcards, output: Path(output[0]).name,
        others=f"--ignore-case --between-stops --max 10000 --min {min(min_probe_size//3,50)}",
    log:
        outdir / "logs/translated_probes/orfipy.log",
    threads: 1
    conda:
        "../../envs/translating.yaml"
    shell:
        """
        orfipy \
        --procs {threads} \
        --pep {params.pep} \
        --outdir {params.out_folder} \
        {params.others} {input} && \ 
        mv {params.out_folder}/*.log {log}
        """


# ToDo: Review the current parsing, with too many orfs it is taking forever
rule parse_translation:
    input:
        pep=outdir / "translated_probes/orfipy/longest_orfs.pep",
    output:
        outdir / "translated_probes/longest_cds.fasta",
    log:
        outdir / "logs/translated_probes/trans_dec_parser.log",
    threads: 1
    run:
        transdecoder_parser(str(input.pep), str(output))
