rule translate_probes:
    input:
        probes_path.resolve(),
    output:
        "<results>/longest_cds.fasta",
    params:
        translated_prop=lookup("probes_translation/translated_prop", within=pipeline),
        stop_per_1Kbp=lookup("probes_translation/stop_per_1Kbp", within=pipeline),
        min_exon_length=lookup("probes_translation/min_exon_length", within=pipeline),
    log:
        "<logs>/<stage>/translated_cds.log",
    threads: 1
    conda:
        "../../envs/translating.yaml"
    script:
        "../../../src/run_orfipy.py"
