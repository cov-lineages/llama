"""
"snakemake --nolock --snakefile {input.snakefile:q} "
                        "{params.force} "
                        "{params.quiet_mode} "
                        "--directory {params.tempdir:q} "
                        "--config "
                        "tempdir={params.tempdir:q} "
                        "reference_fasta={input.reference_fasta:q} "
                        "seqs={input.seqs:q} "
                        "metadata={input.metadata:q} "
                        "to_find_closest={output.to_find_closest:q} "
                        "search_field={params.search_field} "
                        "trim_start={params.trim_start} "
                        "trim_end={params.trim_end} "
                        "--cores {params.cores}"
"""
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
    shell:
        """
        minimap2 -a -x asm5 {input.reference:q} {input.fasta:q} > {output.sam:q}
        """

rule datafunk_trim_and_pad:
    input:
        sam = rules.minimap2_to_reference.output.sam,
        reference = config["reference_fasta"]
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
    output:
        paf = os.path.join(config["tempdir"],"post_qc_query.mapped.paf")
    shell:
        """
        minimap2 -x asm5  -t {threads} --secondary=no --paf-no-hit {input.seqs:q} {input.query_seqs:q} > {output.paf:q}
        """

rule parse_paf:
    input:
        paf = rules.minimap2_against_all.output.paf,
        metadata = config["metadata"],
        fasta = config["seqs"]
    params:
        search_field = config["search_field"]
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
        --search-field {params.search_field}
        """
