"""
Passed into config:

catchment_str=tree_1,tree_2,...,tree_X
"snakemake --nolock --snakefile {input.snakefile_collapse_before:q} "
                            "{params.force} "
                            "{params.quiet_mode} "
                            # "--directory {params.tempdir:q} "
                            "--config "
                            f"catchment_str={catchment_str} "
                            "outdir={params.outdir:q} "
                            # "tempdir={params.tempdir:q} "
                            "not_cog_csv={input.not_cog_csv:q} "
                            "post_qc_query={input.not_cog_query_seqs:q} "
                            "all_cog_seqs={input.all_cog_seqs:q} "
                            "combined_metadata={input.combined_metadata:q} "
                            "--cores {params.cores}"
"""
from Bio import Phylo
from Bio import SeqIO
import csv
import collections

config["tree_stems"] = config["local_str"].split(",")

rule all:
    input:
        expand(os.path.join(config["tempdir"], "collapsed_trees","{tree}.newick"), tree = config["tree_stems"]),
        os.path.join(config["outdir"],"local_trees","collapse_report.txt"),
        expand(os.path.join(config["outdir"],"local_trees","{tree}.tree"), tree = config["tree_stems"])

rule extract_taxa_from_catchment:
    input:
        catchment_tree = os.path.join(config["outdir"],"catchment_trees","{tree}.newick")
    output:
        tree_taxa = os.path.join(config["tempdir"], "catchment_trees","{tree}_taxon_names.txt")
    run:
        print(f"Extracting tax labels from local tree {input.catchment_tree}\n")
        shell("clusterfunk get_taxa -i {input.catchment_tree} --in-format newick -o {output.tree_taxa} --out-format newick")

rule get_lineage_represenatives:
    input:
        metadata = config["metadata"],
        seqs = config["seqs"],
        tree_taxa = rules.extract_taxa_from_catchment.output.tree_taxa
    params:
        data_column = config["data_column"]
    output:
        representative_metadata = os.path.join(config["tempdir"], "representative_lineage_taxa","{tree}.metadata.csv"),
    run:
        taxa = []
        with open(input.tree_taxa, "r") as f:
            for l in f:
                l = l.rstrip("\n")
                taxa.append(l)
        print(f"{len(taxa)} taxa read in")

        lineages = collections.defaultdict(list)
        with open(input.metadata,newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row[params.data_column] in taxa:
                    lineage = row["lineage"]
                    lineages[lineage].append(row[params.data_column])
        print(f"{len(lineages)} different lineages in local tree")
        for lineage in sorted(lineages):
            print(f"- {lineage}")

        lineage_seqs_with_ambiguities = collections.defaultdict(list)
        for record in SeqIO.parse(input.seqs,"fasta"):
            for lineage in lineages:
                if record.id in lineages[lineage]:
                    amb_count = 0
                    for base in record.seq:
                        if base.upper() not in ["A","T","C","G","-"]:
                            amb_count +=1
                    amb_pcent = (100*amb_count) / len(record.seq)

                    lineage_seqs_with_ambiguities[lineage].append((record.id, amb_pcent))

        with open(output.representative_metadata, "w") as fw:
            for lineage in lineage_seqs_with_ambiguities:
                records = lineage_seqs_with_ambiguities[lineage]
                sorted_with_amb = sorted(records, key = lambda x : x[1])
                
                if len(sorted_with_amb) > 10:
                    top_ten_rep = sorted_with_amb[:10]
                    for rep in top_ten_rep:
                        fw.write(f"{rep[0]},{lineage}\n")
                else:
                    for rep in sorted_with_amb:
                        fw.write(f"{rep[0]},{lineage}\n")

rule combine_protected_metadata:
    input:
        lineage_reps = rules.get_lineage_represenatives.output.representative_metadata,
        query_metadata = config["combined_metadata"]
    output:
        protected = os.path.join(config["tempdir"],"collapsed_trees","{tree}.protected.csv")
    run:
        with open(output.protected, "w") as fw:
            fw.write("taxon,lineage\n")
            taxa = []
            with open(input.query_metadata,newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    closest = row["closest"]
                    lineage = row["lineage"]
                    taxa.append(closest)
                    fw.write(f"{closest},{lineage}\n")

            with open(input.lineage_reps,"r") as f:
                for l in f:
                    l = l.rstrip("\n")
                    taxon,lineage = l.split(',')
                    if taxon not in taxa:
                        fw.write(f"{l}\n")

rule summarise_polytomies:
    input:
        tree = os.path.join(config["outdir"], "catchment_trees","{tree}.newick"),
        metadata = rules.combine_protected_metadata.output.protected
    params:
        tree_dir = os.path.join(config["outdir"],"catchment_trees"),
        threshold = config["threshold"]
    output:
        collapsed_tree = os.path.join(config["tempdir"],"collapsed_trees","{tree}.newick"),
        collapsed_information = os.path.join(config["outdir"],"local_trees","{tree}.txt")
    shell:
        """
        clusterfunk focus -i {input.tree:q} \
        -o {output.collapsed_tree:q} \
        --metadata {input.metadata:q} \
        --index-column taxon \
        --in-format newick \
        --out-format newick \
        --threshold {params.threshold} \
        --output-tsv {output.collapsed_information:q}
        """

rule remove_str_for_baltic:
    input:
        tree = os.path.join(config["tempdir"],"collapsed_trees","{tree}.newick")
    output:
        tree = os.path.join(config["outdir"],"local_trees","{tree}.tree")
    run:
        with open(output.tree,"w") as fw:
            with open(input.tree, "r") as f:
                for l in f:
                    l = l.rstrip("\n")
                    l = l.replace("'","")
                    fw.write(l)

rule summarise_processing:
    input:
        collapse_reports = expand(os.path.join(config["outdir"],"local_trees","{tree}.txt"), tree=config["tree_stems"])
    output:
        report = os.path.join(config["outdir"],"local_trees","collapse_report.txt")
    run:
        with open(output.report, "w") as fw:
            for report in input.collapse_reports:
                fn = os.path.basename(report)
                with open(report, "r") as f:
                    for l in f:
                        l = l.rstrip("\n")
                        new_l = f"{fn}\t{l}\n"
                        fw.write(new_l)
