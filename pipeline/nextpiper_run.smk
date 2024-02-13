SAMPLES, = glob_wildcards("reads/pe/{sample}1.fq.gz")
GROUPS = ["1", "2"]
# Find all conda env directories: conda config --show envs_dirs
shell("conda config --add envs_dirs /home/yjkbertrand/miniforge3/envs")
#print(SAMPLES)
#import os
#SAMPLES = ["H1_A1_R", "H1_A2_R"]
#from pathlib import Path
#import os
#configfile: "config.yaml"

#SAMPLE, = glob_wildcards("reads/pe/{sample}1.fq.gz")
#SAMPLE = [x.split("1.fq.gz")[0] for x in os.listdir("./reads/pe/") if x.endswith("1.fq.gz")]

rule all:
    input:
        #trimmed=expand(["trimmed/pe/{sample}1.fq.gz", "trimmed/pe/{sample}2.fq.gz"], sample=SAMPLES)
        #expand("report/pe/{sample}.html", sample=SAMPLES)
        expand("qc/fastqc/{sample}{group}.html", sample=SAMPLES, group=GROUPS)
        #expand("hybpiper_assembled/{sample}/", sample=SAMPLES)

# Read trimming, quality filtering and adaptor removal
rule fastp_pe:
    input:
        sample=["reads/pe/{sample}1.fq.gz", "reads/pe/{sample}2.fq.gz"]
    output:
        trimmed=["trimmed/pe/{sample}1.fq.gz", "trimmed/pe/{sample}2.fq.gz"],
        unpaired=temporary("trimmed/pe/{sample}.singletons.fastq"),
        #failed="trimmed/pe/{sample}.failed.fastq",
        html="report/pe/{sample}.html",
        json=temporary("report/pe/{sample}.json")
    log:
        "logs/fastp/pe/{sample}.log"
    threads: 4
    params:
        extra="--dedup"

    wrapper:
        "v2.2.1/bio/fastp"

# Getting read statistics
rule fastqc:
    input:
        ["trimmed/pe/{sample}{group}.fq.gz"]
        #expand("trimmed/pe/reads/{sample}{group}.fastq", sample=SAMPLES, group = GROUPS)
    output:
        html="qc/fastqc/{sample}{group}.html",
        zip="qc/fastqc/{sample}{group}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}{group}.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v2.2.1/bio/fastqc"

# Hybpiper assembly phase
rule hybpiper_assemble:
    input:
        trimmed_fq = ["trimmed/pe/{sample}1.fq.gz", "trimmed/pe/{sample}2.fq.gz"]
    output:
        assemble=directory("hybpiper_assembled/")
    threads: 8
    log:
        temporary("hybpiper_{sample}.log")
    conda:
      "hybpiper"
    shell:
        "hybpiper assemble -t_dna target_file/targets.fasta -r {input.trimmed_fq} --prefix {wildcards.sample} --bwa --cpu {threads} --keep_intermediate_files --hybpiper_output {output.assemble}"


rule hybpiper_stats:
#hybpiper stats
# Alignment building phase
