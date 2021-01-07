#!/usr/bin/env python3
import csv
import argparse
from Bio import SeqIO
import collections
"""
--csv {input.csv:q} 
--metadata {input.metadata:q} 
--csv-out {output.csv:q} 
"""
def parse_args():
    parser = argparse.ArgumentParser(description='Parse barcode info and csv file, create report.')

    parser.add_argument("--csv", action="store", type=str, dest="input_csv")
    parser.add_argument("--metadata", action="store", type=str, dest="metadata")
    parser.add_argument("--data-column", action="store",type=str, dest="data_column")
    parser.add_argument("--csv-out", action="store", type=str, dest="outfile")
    parser.add_argument("--seqs", action="store", type=str, dest="seqs")
    parser.add_argument("--seqs-out", action="store", type=str, dest="seqs_out")
    return parser.parse_args()


def get_closest_sequences(input_csv):

    closest_cog_sequences = []
    closest_to_query = collections.defaultdict(list)
    with open(input_csv,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            closest_to_query[row["closest"]].append(row)

    return closest_to_query


def parse_csv_and_get_metadata():
    args = parse_args()

    closest_to_query = get_closest_sequences(args.input_csv)
    column_to_match = args.data_column
    
    with open(args.metadata, newline="") as f:
        rows_to_write = []
        reader = csv.DictReader(f)
        header_names = reader.fieldnames
        with open(args.outfile, "w") as fw:
            header_names.append("query")
            header_names.append("closest")
            header_names.append("SNPdistance")
            header_names.append("SNPs")
            writer = csv.DictWriter(fw, fieldnames=header_names,lineterminator='\n')
            writer.writeheader()
        
            for row in reader:
                if row[column_to_match] in closest_to_query:
                    for query_row in closest_to_query[row[column_to_match]]:
                        new_row = row
                        new_row["query"]=query_row["query"]
                        new_row["closest"]=row[column_to_match]
                        new_row["SNPdistance"]=query_row["SNPdistance"]
                        new_row["SNPs"]=query_row["SNPs"]

                        writer.writerow(new_row)


    with open(args.seqs_out, "w") as fw:
        for record in SeqIO.parse(args.seqs,"fasta"):
            if record.id in closest_to_query:
                queries = []
                for query_row in closest_to_query[record.id]:
                    queries.append(query_row["query"])
                closest_queries = ",".join(queries)
                fw.write(f">{record.id} query={closest_queries}\n{record.seq}\n")

if __name__ == '__main__':

    parse_csv_and_get_metadata()
    