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
if config["report"]==True:
    rule all:
        input:
            os.path.join(config["outdir"],"catchment_trees","catchment_trees_prompt.txt"),
            os.path.join(config["outdir"],"not_in_db.csv"),
            os.path.join(config["outdir"],"local_trees","collapse_report.txt"),
            os.path.join(config["outdir"], "llama_report.md")
else:
    rule all:
        input:
            os.path.join(config["outdir"],"catchment_trees","catchment_trees_prompt.txt"),
            os.path.join(config["outdir"],"not_in_db.csv"),
            os.path.join(config["outdir"],"local_trees","collapse_report.txt")


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
        txt = os.path.join(config["outdir"],"catchment_trees","catchment_trees_prompt.txt")
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

rule process_local_trees:
    input:
        snakefile_just_collapse = os.path.join(workflow.current_basedir,"just_collapse_no_seqs.smk"), #alternative snakefiles
        combined_metadata = rules.check_metadata.output.in_db, 
        context_prompt = rules.jclusterfunk_context.output.txt,
        metadata = config["metadata"]
    params:
        outdir= config["outdir"],
        tempdir= config["tempdir"],
        path = workflow.current_basedir,
        threshold = config["threshold"],
        data_column = config["data_column"],
        lineage_representatives = config["lineage_representatives"],
        number_of_representatives = config["number_of_representatives"],
        tree_dir = os.path.join(config["outdir"],"catchment_trees"),
        cores = workflow.cores,
        force = config["force"],
        quiet_mode = config["quiet_mode"]
    output:
        tree_summary = os.path.join(config["outdir"],"local_trees","collapse_report.txt")
    run:
        local_trees = []
        for r,d,f in os.walk(params.tree_dir):
            for fn in f:
                if fn.endswith(".newick"):
                    file_stem = ".".join(fn.split(".")[:-1])
                    local_trees.append(file_stem)
        local_str = ",".join(local_trees) #to pass to snakemake pipeline

        shell("snakemake --nolock --snakefile {input.snakefile_just_collapse:q} "
                        "{params.force} "
                        "{params.quiet_mode} "
                        "--directory {params.tempdir:q} "
                        "--config "
                        f"local_str={local_str} "
                        "outdir={params.outdir:q} "
                        "tempdir={params.tempdir:q} "
                        "number_of_representatives={params.number_of_representatives} "
                        "lineage_representatives={params.lineage_representatives} "
                        "metadata={input.metadata:q} "
                        "data_column={params.data_column} "
                        "combined_metadata={input.combined_metadata:q} "
                        "threshold={params.threshold} "
                        "--cores {params.cores}")


rule make_report:
    input:
        lineage_trees = rules.process_local_trees.output.tree_summary,
        query = config["query"],
        combined_metadata = rules.check_metadata.output.in_db,
        metadata = config["metadata"],
        footer = config["footer"],
        report_template = config["report_template"],
        no_seq = rules.check_metadata.output.not_in_db
    params:
        treedir = os.path.join(config["outdir"],"local_trees"),
        outdir = config["rel_outdir"],
        rel_figdir = os.path.join(".","figures"),
        figdir = os.path.join(config["outdir"],"figures"),
        failure = config["qc_fail"],
        input_column = config["input_column"],
        data_column = config["data_column"],
        colour_fields = config["colour_fields"],
        label_fields = config["label_fields"]
    output:
        outfile = os.path.join(config["outdir"], "llama_report.md"), 
        footer_fig = os.path.join(config["outdir"], "figures", "footer.png")
    run:
        shell("""
            cp {input.footer:q} {output.footer_fig:q}
        """)
        shell(
            "make_report.py "
            "--input-csv {input.query:q} "
            "--figdir {params.rel_figdir:q} "
            "{params.failure} "
            "--no-seq-provided {input.no_seq} "
            "--treedir {params.treedir:q} "
            "--report-template {input.report_template:q} "
            "--filtered-metadata {input.combined_metadata:q} "
            "--metadata {input.metadata:q} "
            "--outfile {output.outfile:q} "
            "--outdir {params.outdir:q} "
            "--input-column {params.input_column:q} "
            "--data-column {params.data_column:q} "
            "--colour-fields {params.colour_fields:q} "
            "--label-fields {params.label_fields:q}")
