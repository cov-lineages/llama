import csv
from Bio import SeqIO
import os
import collections

rule check_metadata:
    input:
        query = config["query"],
        seqs = config["seqs"],
        metadata = config["metadata"]
    message: "Checking metadata for queries"
    params:
        input_column = config["input_column"],
        data_column = config["data_column"]
    output:
        in_db = os.path.join(config["tempdir"],"query_in_db.csv"),
        seqs = os.path.join(config["tempdir"],"query_in_db.fasta"),
        not_in_db = os.path.join(config["tempdir"],"not_in_db.csv")
    shell:
        """
        check_metadata.py --query {input.query:q} \
                        --seqs {input.seqs:q} \
                        --metadata {input.metadata:q} \
                        --input-column {params.input_column} \
                        --data-column {params.data_column} \
                        --in-metadata {output.in_db:q} \
                        --in-seqs {output.seqs:q} \
                        --not-in-db {output.not_in_db:q}
        """

rule get_closest_in_db:
    input:
        snakefile = os.path.join(workflow.current_basedir,"find_closest_in_db.smk"),
        reference_fasta = config["reference_fasta"],
        seqs = config["seqs"],
        metadata = config["metadata"],
        not_in_db = rules.check_metadata.output.not_in_db, #use
    params:
        outdir= config["outdir"],
        tempdir= config["tempdir"],
        path = workflow.current_basedir,
        cores = workflow.cores,
        force = config["force"],
        fasta = config["fasta"], #use
        input_column = config["input_column"],
        data_column = config["data_column"],
        query = config["post_qc_query"], #use
        stand_in_query = os.path.join(config["tempdir"], "temp.fasta"),
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
        quiet_mode = config["quiet_mode"]
    output:
        closest_in_db = os.path.join(config["tempdir"],"closest_in_db.csv"),
        to_find_closest = os.path.join(config["tempdir"],"to_find_closest.fasta"),
        aligned_query = os.path.join(config["tempdir"],"post_qc_query.aligned.fasta"),
        not_processed = os.path.join(config["tempdir"], "no_seq_to_process.csv")
    run:
        query_with_no_seq = []
        to_find_closest = {}

        not_db = []
        with open(input.not_in_db, newline = "") as f: # getting list of non-db queries
            reader = csv.DictReader(f)
            for row in reader:
                not_db.append(row[params.input_column])

        if params.fasta != "":
             # get set with supplied sequences
                print("Not in db but have a sequence supplied:")
                for record in SeqIO.parse(params.query, "fasta"):
                    if record.id in not_db:
                        print(f"\t-{record.id}")
                        to_find_closest[record.id] = record.seq # overwrites with supplied seq if found in all db

        with open(output.to_find_closest, "w") as fw:
            for seq in to_find_closest:
                fw.write(f">{seq}\n{to_find_closest[seq]}\n")

        for query in not_db: # get set with no sequences supplied
            if query not in to_find_closest:
                query_with_no_seq.append(query)

        if to_find_closest != {}:
            print(f"Passing {len(to_find_closest)} sequences into find_closest_in_db search pipeline")
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                        "{params.force} "
                        "{params.quiet_mode} "
                        "--directory {params.tempdir:q} "
                        "--config "
                        "tempdir={params.tempdir:q} "
                        "reference_fasta={input.reference_fasta:q} "
                        "seqs={input.seqs:q} "
                        "metadata={input.metadata:q} "
                        "to_find_closest={output.to_find_closest:q} "
                        "data_column={params.data_column} "
                        "trim_start={params.trim_start} "
                        "trim_end={params.trim_end} "
                        "--cores {params.cores}")

        else:
            shell("touch {output.closest_in_db:q} && touch {output.aligned_query} ")
        print(f"Number of query ids not processed: {len(query_with_no_seq)}\n Please provide a valid fasta sequence if you'd like these to be included.")
        with open(output.not_processed, "w") as fw:
            for query in list(set(query_with_no_seq)):
                fw.write(f"{query},fail=no_sequence\n")

rule combine_metadata:
    input:
        closest_in_db = rules.get_closest_in_db.output.closest_in_db,
        in_db = rules.check_metadata.output.in_db
    params:
        data_column = config["data_column"]
    message: "Combining metadata from those already in db and the closest ones in db"
    output:
        combined_csv = os.path.join(config["outdir"],"combined_metadata.csv")
    run:
        with open(output.combined_csv,"w") as fw:
            with open(input.in_db, "r") as f:
                for l in f:
                    l = l.rstrip("\n")
                    fw.write(l + '\n')
            with open(input.closest_in_db, "r") as f:
                for l in f:
                    l = l.rstrip("\n")
                    if params.data_column in l:
                        pass
                    else:
                        fw.write(l + '\n')

rule jclusterfunk_context:
    input:
        tree = config["tree"],
        metadata = rules.combine_metadata.output.combined_csv
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
        snakefile = os.path.join(workflow.current_basedir,"process_local_trees.smk"), #alternative snakefiles
        snakefile_just_collapse = os.path.join(workflow.current_basedir,"just_collapse_trees.smk"), #alternative snakefiles
        combined_metadata = rules.combine_metadata.output.combined_csv, 
        query_seqs = rules.get_closest_in_db.output.aligned_query, #datafunk-processed seqs
        context_prompt = rules.jclusterfunk_context.output.txt,
        metadata = config["metadata"],
        seqs = config["seqs"],
        reference_fasta = config["reference_fasta"]
    params:
        outdir= config["outdir"],
        tempdir= config["tempdir"],
        path = workflow.current_basedir,
        threshold = config["threshold"],
        data_column = config["data_column"],
        fasta = config["fasta"],
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

        query_seqs = 0
        for record in SeqIO.parse(input.query_seqs,"fasta"):
            query_seqs +=1

        if query_seqs !=0:
            print(f"Passing {input.query_seqs} into processing pipeline.")
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                        "{params.force} "
                        "{params.quiet_mode} "
                        "--directory {params.tempdir:q} "
                        "--config "
                        f"local_str={local_str} "
                        "outdir={params.outdir:q} "
                        "tempdir={params.tempdir:q} "
                        "outgroup_fasta={input.reference_fasta} "
                        "aligned_query_seqs={input.query_seqs:q} "
                        "metadata={input.metadata:q} "
                        "seqs={input.seqs:q} "
                        "data_column={params.data_column} "
                        "combined_metadata={input.combined_metadata:q} "
                        "threshold={params.threshold} "
                        "--cores {params.cores}")
        else:
            print(f"No new sequences to add in, just collapsing trees.")
            shell("snakemake --nolock --snakefile {input.snakefile_just_collapse:q} "
                            "{params.force} "
                            "{params.quiet_mode} "
                            "--directory {params.tempdir:q} "
                            "--config "
                            f"local_str={local_str} "
                            "outdir={params.outdir:q} "
                            "tempdir={params.tempdir:q} "
                            "outgroup_fasta={input.reference_fasta} "
                            "aligned_query_seqs={input.query_seqs:q} "
                            "metadata={input.metadata:q} "
                            "seqs={input.seqs:q} "
                            "data_column={params.data_column} "
                            "combined_metadata={input.combined_metadata:q} "
                            "threshold={params.threshold} "
                            "--cores {params.cores}")
