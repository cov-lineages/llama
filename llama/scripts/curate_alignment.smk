import csv
from Bio import SeqIO
import os
import collections



rule minimap2_db_to_reference:
    input:
        fasta = config["seqs"],
        reference = config["reference_fasta"]
    output:
        sam = os.path.join(config["tempdir"],"post_qc_aln.reference_mapped.sam")
    message: "Running minimap2 against the reference (early lineage A) sequence"
    log:
        os.path.join(config["tempdir"],"logs","minimap2_to_reference.log")
    shell:
        """
        minimap2 -a -x asm5 {input.reference:q} {input.fasta:q} -o {output.sam:q} &> {log}
        """

rule datafunk_trim_and_pad:
    input:
        sam = rules.minimap2_db_to_reference.output.sam,
        reference = config["reference_fasta"]
    message: "Running datafunk to trim and pad against the reference"
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
        insertions = os.path.join(config["tempdir"],"post_qc_query.insertions.txt")
    output:
        fasta = os.path.join(config["tempdir"],"seqs_db.aligned.fasta")
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
