#!/usr/bin/env python

import csv
from Bio import SeqIO
import os
import collections

##### Configuration #####

if config.get("force"):
    config["force"] = "--forceall "

if config["quiet_mode"]=="True":
    config["quiet_mode"] = "--quiet"
else:
    config["quiet_mode"] = ""

config["distance"]=int(config["distance"])

##### Target rules #####
if config["report"]=="True":
    rule all:
        input:
            os.path.join(config["tempdir"],"catchment_trees","catchment_trees_prompt.txt"),
            os.path.join(config["outdir"],"not_in_db.csv")
else:
    rule all:
        input:
            os.path.join(config["tempdir"],"catchment_trees","catchment_trees_prompt.txt"),
            os.path.join(config["outdir"],"not_in_db.csv")


rule check_metadata:
    input:
        query = config["query"],
        metadata = config["metadata"]
    message: "Checking metadata for queries"
    params:
        input_column = config["input_column"],
        data_column = config["data_column"]
    output:
        in_db = os.path.join(config["tempdir"],"query_in_db.csv"),
        not_in_db = os.path.join(config["outdir"],"not_in_db.csv")
    run:
        
        found = []
        in_metadata = []
        in_metadata_names = []
        
        query_names = []
        with open(input.query,newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                query_names.append(row[params.input_column])

        with open(input.metadata,newline="") as f:
            reader = csv.DictReader(f)
            header_names = reader.fieldnames
            for row in reader:
                for seq in query_names:

                    db_id = row[params.data_column]
                    if seq == db_id:
                        
                        row["query"]=row[params.data_column]
                        row["closest"]=row[params.data_column]
                        in_metadata.append(row)
                        in_metadata_names.append(seq)

        print(f"Number of query records found in metadata: {len(in_metadata)}")
        with open(output.in_db, "w") as fw:
            header_names.append("query")
            header_names.append("closest")
            writer = csv.DictWriter(fw, fieldnames=header_names,lineterminator='\n')
            writer.writeheader()
            writer.writerows(in_metadata)

        print(f"Number of associated sequences found: {len(found)}")
        with open(output.not_in_db, "w") as fw:
            print("\nThe following sequences were not found in the database:")
            fw.write(f"{params.input_column}\n")
            for query in query_names:
                if query not in found:
                    fw.write(query + '\n')
                    print(f"\t-{query}")

rule jclusterfunk_context:
    input:
        tree = config["tree"],
        metadata = rules.check_metadata.output.in_db
    message: "Finding the local context of query sequences in the global tree"
    params:
        outdir = os.path.join(config["outdir"],"catchment_trees"),
        distance = config["distance"]
    output:
        txt = os.path.join(config["tempdir"],"catchment_trees","catchment_trees_prompt.txt")
    shell:
        """
        jclusterfunk context \
        -i {input.tree:q} \
        -o {params.outdir:q} \
        --max-parent {params.distance} \
        -f newick \
        -p local \
        --ignore-missing \
        -m {input.metadata:q} \
        --id-column closest \
        && touch {output.txt:q} 
        """

