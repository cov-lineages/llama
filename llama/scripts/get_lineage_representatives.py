#!/usr/bin/env python3

import argparse
import collections
from Bio import SeqIO
import os
import csv
cwd = os.getcwd()


def parse_args():
    parser = argparse.ArgumentParser(description='Get selection of lineage representatives.')

    parser.add_argument("--tree-taxa", action="store", type=str, dest="tree_taxa")
    parser.add_argument("--seqs", action="store", type=str, dest="seqs")
    parser.add_argument("--metadata", action="store", type=str, dest="metadata")
    parser.add_argument("--data-column", action="store", type=str, dest="data_column")
    parser.add_argument("--representatives", action="store", type=str, dest="representatives")
    parser.add_argument("--number-of-representatives", action="store", type=int, dest="number", default=5)
    return parser.parse_args()

def get_representatives():
    args = parse_args()

    taxa = []
    with open(args.tree_taxa, "r") as f:
        for l in f:
            l = l.rstrip("\n")
            taxa.append(l)
    print(f"{len(taxa)} taxa read in")

    lineages = collections.defaultdict(list)
    with open(args.metadata,newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row[args.data_column] in taxa:
                lineage = row["lineage"]
                lineages[lineage].append(row[args.data_column])
    print(f"{len(lineages)} different lineages in local tree")
    for lineage in sorted(lineages):
        print(f"- {lineage}")

    lineage_seqs_with_ambiguities = collections.defaultdict(list)
    for record in SeqIO.parse(args.seqs,"fasta"):
        for lineage in lineages:
            if record.id in lineages[lineage]:
                amb_count = 0
                for base in record.seq:
                    if base.upper() not in ["A","T","C","G","-"]:
                        amb_count +=1
                amb_pcent = (100*amb_count) / len(record.seq)

                lineage_seqs_with_ambiguities[lineage].append((record.id, amb_pcent))

    num = int(args.number)
    with open(args.representative_metadata, "w") as fw:
        for lineage in lineage_seqs_with_ambiguities:
            records = lineage_seqs_with_ambiguities[lineage]
            sorted_with_amb = sorted(records, key = lambda x : x[1])
            

            if len(sorted_with_amb) > num:
                top_ten_rep = sorted_with_amb[:num]
                for rep in top_ten_rep:
                    fw.write(f"{rep[0]},{lineage}\n")
            else:
                for rep in sorted_with_amb:
                    fw.write(f"{rep[0]},{lineage}\n")

if __name__ == '__main__':

    get_representatives()