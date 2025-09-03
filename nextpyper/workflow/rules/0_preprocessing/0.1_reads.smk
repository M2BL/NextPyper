rule fastp_pe:
    input:
        in1=lookup(query="sample=='{sample}'", cols="path_forward", within=sample_table),
        in2=lookup(query="sample=='{sample}'", cols="path_reverse", within=sample_table),
    output:
        trim1=outdir / "preprocessed/trimmed/{sample}_R1.fastq.gz",
        trim2=outdir / "preprocessed/trimmed/{sample}_R2.fastq.gz",
        html=outdir / "logs/preprocessing/fastp/{sample}.html",
        json=outdir / "logs/preprocessing/fastp/{sample}.json",
    log:
        outdir / "logs/preprocessing/fastp/{sample}.log",
    params:
        trim_qual=lookup("preprocessing/trim_qual", within=pipeline),
        trim_min_len=lookup("preprocessing/min_len", within=pipeline),
        extra="--trim_poly_g --trim_poly_x --low_complexity_filter --cut_tail",
    threads: 4
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """fastp --thread {threads} {params.extra} \
        --cut_mean_quality {params.trim_qual} \
        --length_required {params.trim_min_len} \
        --in1 {input.in1} --in2 {input.in2} \
        --out1 {output.trim1} --out2 {output.trim2} \
        --html {output.html} \
        --json {output.json} 2> {log} 
        """


checkpoint prepare_cps:
    output:
        outdir / "preprocessed/ref_cps.fasta",
    log:
        outdir / "logs/preprocessing/cps.log",
    params:
        cp_refs_map=cp_refs_map,
        custom=lookup("args/custom_cps", within=config),
    retries: 5
    run:
        out_cps = Path(output[0])
        custom_cps = params.custom

        with open(log[0], "w") as outlog:
            sys.stdout = sys.stderr = outlog
            cps = []

            if use_ref_cps and seqids:
                handle = Entrez.efetch(
                    db="nuccore", id=seqids, rettype="fasta", retmode="text"
                )
                cps += list(SeqIO.parse(handle, "fasta"))

            # If custom cps are given add them to the downloaded ones
            if custom_cps:
                cps += list(SeqIO.parse(custom_cps, "fasta"))

            SeqIO.write(cps, out_cps, "fasta")


rule reads_cp_cleaning:
    input:
        ref=outdir / "preprocessed/ref_cps.fasta",
        trim1=outdir / "preprocessed/trimmed/{sample}_R1.fastq.gz",
        trim2=outdir / "preprocessed/trimmed/{sample}_R2.fastq.gz",
    output:
        clean1=outdir / "preprocessed/cleaned/{sample}_R1.fastq.gz",
        clean2=outdir / "preprocessed/cleaned/{sample}_R2.fastq.gz",
    log:
        outdir / "logs/preprocessing/cleaning/{sample}.log",
    conda:
        "../../envs/preprocessing.yaml"
    threads: 4
    run:
        if Path(input.ref).stat().st_size > 0:
            shell(
                """minimap2 -t {threads} -ax sr {input.ref} {input.trim1} {input.trim2} 2> {log} | \
                        samtools view -u -e 'flag & 4 || flag & 8' 2> {log} | \
                        samtools fastq -1 {output.clean1} -2 {output.clean2} 2> {log} """
            )
        else:
            Path(output.clean1).symlink_to(
                Path(input.trim1).relative_to(Path(output.clean1).parent, walk_up=True)
            )
            Path(output.clean2).symlink_to(
                Path(input.trim2).relative_to(Path(output.clean2).parent, walk_up=True)
            )
