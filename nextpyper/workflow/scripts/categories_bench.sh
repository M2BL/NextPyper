#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

## DESCRIPTION:
## This script computes the categories benchmark for a run comparing a set of queries 
## (a program output) and a set of targets (a gold standard). The categorization is further
## explained in bench_busco.py.  
## Keeping hits with at least 25% coverage

## REQUIREMENTS:
## Magicblast, Vsearch, GNU parallel, python>=3.10 (with Biopython, polars and intervaltree) 
## should be available in the PATH. bench_busco.py should be in the same directory as this script.
## The script will create a blast database for the targets if it does not exist.

## POSITIONAL COMMAND LINE ARGUMENTS:
## 1. Targets directory - gold standard targets in a single folder. 
## 2. Queries directory - Sequence samples to be evaluated (in fasta format).
## 3. Tables directory - Path to write the output tables
## 4. Max threads to use (the number of concurrent jobs will threads/4)

## EXAMPLE OF USAGE:
## mmseqs_bench.sh <targets_dir> <queries_dir> <tables_dir> <max_threads>

targets_dir=$1
queries_dir=$2
tables_dir=$3
max_threads=$4
jobs=$((max_threads / 4))


function blast_create_db {
    local targets=$1
    local db_name=$2
    makeblastdb -dbtype nucl -parse_seqids -in "${targets}" -out "${db_name}"
}

function magicblast_match {
    ## POSITIONAL COMMAND LINE ARGUMENTS:
    ## 1. Targets db - targets as a blast database. 
    ## 2. Path to assembled sequences (query) - in fasta format.
    ## 3. Path to write the output table
    ## 4. Threads to use
    
    local targets=$1
    local asm=$2
    local out_table=$3
    local threads=$4

    magicblast -splice F -outfmt "tabular" \
        -num_threads "${threads}" \
        -query "${asm}" \
        -db "${targets}" \
        -out "${out_table}" 
}

function vsearch_chimera {
    ## Run reference-based chimera detection with vsearch

    local asm=$1
    local targets=$2
    local out_table=$3
    local threads=$4
    vsearch --db "${targets}" --uchime_ref "$asm" --uchimeout "${out_table}" --threads "${threads}"
}


export -f blast_create_db
export -f magicblast_match
export -f vsearch_chimera

mkdir -p "${tables_dir}/matches"
mkdir -p "${tables_dir}/vsearch"

# Check if we have to create the blast databases for the targets
orig_targets_dir="${targets_dir}"
array=($(find ${targets_dir} -mindepth 1 -maxdepth 1 -type f -name "*.ndb"))
if [ "${#array[@]}" -eq 0 ]; then
    # We need to created the mmeseqs2 dbs
    mkdir -p "${targets_dir}_dbs"

    array=($(find ${targets_dir}_dbs -mindepth 1 -maxdepth 1 -type f -name "*.ndb"))
    if [ "${#array[@]}" -eq 0 ]; then
        # The dbs really do not exist
        # parallel "mmseqs_create_db {} {//}_dbs/{/.}" ::: "${targets_dir}"/*.fasta
        parallel "blast_create_db {} {//}_dbs/{/.}" ::: "${targets_dir}"/*.fasta
    fi

    targets_dir="${targets_dir}_dbs"
fi

parallel --jobs "${jobs}" magicblast_match "${targets_dir}/{/.} {} ${tables_dir}/matches/{/.}.tsv 4" ::: "${queries_dir}"/*.fasta
parallel --jobs "${jobs}" vsearch_chimera {} "${orig_targets_dir}/{/} ${tables_dir}/vsearch/{/.}.tsv 4" ::: "${queries_dir}"/*.fasta
python "$(dirname $0)/bench_busco.py" --batch "${tables_dir}/matches" "${tables_dir}/vsearch" "${orig_targets_dir}" "${tables_dir}/categories.tsv"
