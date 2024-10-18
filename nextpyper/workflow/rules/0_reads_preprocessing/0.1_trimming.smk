# The input function map the sample at hand to its input files (specified in the samples table):
def get_raw_input_fastq_r1(wildcards):
    return sample_dict[wildcards.sample]["path_forward"]


def get_raw_input_fastq_r2(wildcards):
    return sample_dict[wildcards.sample]["path_reverse"]


if multi_probes and use_ref_cps:

    rule cps_download:
        output:
            outdir / "preprocessed/ref_cps.fasta",
        log:
            outdir / "logs/preprocessing/ref_cps.fasta",
        params:
            min_probes_cov=lookup("cp_cleaning/min_sp_probes_cov", within=pipeline),
        retries: 5
        run:
            out_cps = Path(output[0])
            with open(log[0], "w") as outlog:
                sys.stdout = sys.stderr = outlog

                kp2seqid = pd.read_csv(cp_refs_map)
                selected_sps = {
                    sp
                    for sp, count in Counter(
                        probe.id.split("_")[0] for probe in probes
                    ).items()
                    if count >= params.min_probes_cov
                }
                seqids = (
                    kp2seqid[kp2seqid["1kp"].isin(selected_sps)]["seqid"]
                    .unique()
                    .tolist()
                )
                if seqids:
                    handle = Entrez.efetch(
                        db="nuccore", id=seqids, rettype="fasta", retmode="text"
                    )
                    cps = list(SeqIO.parse(handle, "fasta"))
                    SeqIO.write(cps, out_cps, "fasta")
                else:
                    out_cps.touch()


else:

    rule no_cps_download:
        output:
            touch(outdir / "preprocessed/ref_cps.fasta"),


rule fastp_pe:
    input:
        in1=get_raw_input_fastq_r1,
        in2=get_raw_input_fastq_r2,
    output:
        trim1=outdir / "preprocessed/trimmed/{sample}_R1.fastq.gz",
        trim2=outdir / "preprocessed/trimmed/{sample}_R2.fastq.gz",
        html=outdir / "logs/preprocessing/fastp/{sample}.html",
        json=outdir / "logs/preprocessing/fastp/{sample}.json",
    log:
        outdir / "logs/preprocessing/fastp/{sample}.log",
    params:
        extra="--trim_poly_g --trim_poly_x --low_complexity_filter --cut_right",
    threads: 4
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        "(fastp --thread {threads} "
        "{params.extra} "
        "--in1 {input.in1} --in2 {input.in2} "
        "--out1 {output.trim1} --out2 {output.trim2} "
        "--html {output.html} "
        "--json {output.json} ) 2> {log} "


def make_cleaning_reference(wildcards, input) -> str:
    refs = [input.rrna]

    if multi_probes and use_ref_cps:
        refs.append(str((outdir / "preprocessed/ref_cps.fasta").resolve()))

    if custom_cps:
        refs.append(str(custom_cps.resolve()))

    return ",".join(refs)


rule cleaning_data:
    input:
        in1=outdir / "preprocessed/trimmed/{sample}_R1.fastq.gz",
        in2=outdir / "preprocessed/trimmed/{sample}_R2.fastq.gz",
        rrna=silva_db,
        ref_cps=outdir / "preprocessed/ref_cps.fasta",
    output:
        out1=outdir / "preprocessed/cleaned/{sample}_R1.fastq.gz",
        out2=outdir / "preprocessed/cleaned/{sample}_R2.fastq.gz",
    log:
        outdir / "logs/preprocessing/bbduk_cleaning/{sample}.log",
    params:
        refs=make_cleaning_reference,
        others=other_bbduk,
        k=19,
    threads: 4
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        "(bbduk.sh {params.others} "
        "in1={input.in1} "
        "in2={input.in2} "
        "ref={params.refs} "
        "out1={output.out1} "
        "out2={output.out2} "
        "k={params.k} "
        "threads={threads} ) 2> {log} "


rule matching_probes:
    input:
        in1=outdir / "preprocessed/cleaned/{sample}_R1.fastq.gz",
        in2=outdir / "preprocessed/cleaned/{sample}_R2.fastq.gz",
        ref=probes_path.resolve(),
    output:
        outm1=outdir / "preprocessed/filtered/{sample}_R1.fastq",
        outm2=outdir / "preprocessed/filtered/{sample}_R2.fastq",
    log:
        outdir / "logs/preprocessing/bbduk_probe_matching/{sample}.log",
    params:
        others=other_bbduk,
        k=bbduk_k,
    threads: 4
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        "(bbduk.sh {params.others} "
        "in1={input.in1} "
        "in2={input.in2} "
        "ref={input.ref} "
        "outm1={output.outm1} "
        "outm2={output.outm2} "
        "k={params.k} "
        "threads={threads} ) 2> {log} "
