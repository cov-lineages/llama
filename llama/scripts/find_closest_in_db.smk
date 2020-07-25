import os
from Bio import SeqIO
rule all:
    input:
        os.path.join(config["tempdir"],"closest_in_db.csv"),
        os.path.join(config["tempdir"],"post_qc_query.aligned.fasta")


rule minimap2_to_reference:
    input:
        fasta = config["to_find_closest"],
        reference = config["reference_fasta"]
    output:
        sam = os.path.join(config["tempdir"],"post_qc_query.reference_mapped.sam")
    message: "Running minimap2 against the reference (early lineage A) sequence"
    log:
        os.path.join(config["tempdir"],"logs","minimap2_to_reference.log")
    shell:
        """
        minimap2 -a -x asm5 {input.reference:q} {input.fasta:q} -o {output.sam:q} &> {log}
        """

rule datafunk_trim_and_pad:
    input:
        sam = rules.minimap2_to_reference.output.sam,
        reference = config["reference_fasta"]
    message: "Running datafunk to trim and pad against the reference"
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
        insertions = os.path.join(config["tempdir"],"post_qc_query.insertions.txt")
    output:
        fasta = os.path.join(config["tempdir"],"post_qc_query.aligned.fasta")
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam:q} \
          -r {input.reference:q} \
          -o {output.fasta:q} \
          -t [{params.trim_start}:{params.trim_end}] \
          --pad \
          --log-inserts 
        """

rule minimap2_against_all:
    input:
        query_seqs = rules.datafunk_trim_and_pad.output.fasta,
        seqs = config["seqs"]
    threads: workflow.cores
    message: "Running minimap2 against the entire fasta db using {threads} threads"
    output:
        paf = os.path.join(config["tempdir"],"post_qc_query.mapped.paf")
    log:
        os.path.join(config["tempdir"],"logs","minimap2_to_all.log")
    shell:
        """
        minimap2 -x asm5  -t {threads} --secondary=no \
        --paf-no-hit {input.seqs:q} {input.query_seqs:q} \
        -o {output.paf:q} &> {log}
        """

rule parse_paf:
    input:
        paf = rules.minimap2_against_all.output.paf,
        metadata = config["metadata"],
        fasta = config["seqs"]
    params:
        data_column = config["data_column"]
    message: "Parsing out the top hit in the database for each query sequence"
    output:
        fasta = os.path.join(config["tempdir"],"closest_in_db.fasta"),
        csv = os.path.join(config["tempdir"],"closest_in_db.csv")
    shell:
        """
        parse_paf.py \
        --paf {input.paf:q} \
        --metadata {input.metadata:q} \
        --seqs {input.fasta:q} \
        --csv-out {output.csv:q} \
        --seqs-out {output.fasta:q} \
        --data-column {params.data_column}
        """
