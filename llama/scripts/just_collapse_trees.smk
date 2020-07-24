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

# rule annotate:
#     input:
#         tree = os.path.join(config["tempdir"],"catchment_trees","{tree}.nexus"),
#         metadata = config["combined_metadata"]
#     output:
#         tree = os.path.join(config["outdir"],"annotated_trees","{tree}.nexus")
#     shell:
#         """
#         ~/Documents/jclusterfunk/release/jclusterfunk_v0.0.1/jclusterfunk annotate \
#         -i {input.tree:q} \
#         -o {output.tree} \
#         -m {input.metadata:q} \
#         -r \
#         --id-column closest \
#         --tip-attributes lineage \
#         -f nexus
#         """

rule extract_taxa:
    input:
        catchment_tree = os.path.join(config["outdir"],"catchment_trees","{tree}.nexus")
    output:
        tree_taxa = os.path.join(config["tempdir"], "catchment_trees","{tree}_taxon_names.txt")
    shell:
        "clusterfunk get_taxa -i {input.catchment_tree} --in-format nexus -o {output.tree_taxa} --out-format nexus"

rule get_lineage_represenatives:
    input:
        metadata = config["metadata"],
        seqs = config["seqs"],
        taxa = rules.extract_taxa.output.tree_taxa
    params:
        data_column = config["data_column"]
    output:
        representative_seq = os.path.join(config["tempdir"], "representative_lineage_taxa","{tree}.fasta"),
        representative_metadata = os.path.join(config["tempdir"], "representative_lineage_taxa","{tree}.metadata.csv"),
    run:
        taxa = []
        with open(input.tree_taxa, "r") as f:
            for l in f:
                l = l.rstrip("\n")
                taxa.append(l)

        lineages = collections.defaultdict(list)
        with open(input.metadata,newline="") as f:
            reader = csv.DictReader(f)
            if row[params.data_column] in taxa:
                lineage = row["lineage"]
                lineages[lineage].append(row[params.data_column])

        lineage_seqs = collections.defaultdict(list)
        for record in SeqIO.parse(input.seqs,"fasta"):
            for lineage in lineages:
                if record.id in lineages[lineage]:
                    lineage_seqs[lineage].append(record)

        with open(output.representative_metadata, "w") as fcsv:
            with open(output.representative_seq, "w") as fw:
                for lineage in lineage_seqs:
                    records = lineage_seqs[lineage]
                    sorted_with_amb = []
                    for record in records:
                        amb_count = 0
                        for base in record.seq:
                            if base.upper() not in ["A","T","C","G","-"]:
                                amb_count +=1
                        amb_pcent = (100*amb_count) / len(record.seq)
                        sorted_with_amb.append((record.id, amb_pcent, record.seq))
                    sorted_with_amb = sorted(sorted_with_amb, key = lambda x : x[1])
                    top_five_rep = sorted_with_amb[:5]
                    for rep in top_five_rep:
                        fw.write(f">{rep[0]} lineage={lineage} ambiguity={rep[1]}\n{rep[2]}\n")
                        fcsv.write(f"{rep[0]},{lineage}\n")

rule combine_protected_metadata:
    input:
        lineage_reps = rules.get_lineage_represenatives.output.representative_metadata,
        query_metadata = config["combined_metadata"]
    output:
        protected = os.path.join(config["tempdir"],"collapsed_trees","{tree}.protected.csv")
    run:
        with open(output.protected, "w") as fw:
            fw.write("taxon,lineage\n")
            with open(input.query_metadata,newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    closest = row["closest"]
                    lineage = row["lineage"]
                    fw.write(f"{closest},{lineage}\n")
            with open(input.lineage_reps,"r") as f:
                for l in f:
                    l = l.rstrip("\n")
                    fw.write(f"{l}\n")

rule summarise_polytomies:
    input:
        tree = os.path.join(config["outdir"], "catchment_trees","{tree}.nexus"),
        metadata = config["combined_metadata"]
    params:
        tree_dir = os.path.join(config["outdir"],"catchment_trees")
    output:
        collapsed_tree = os.path.join(config["tempdir"],"collapsed_trees","{tree}.newick"),
        collapsed_information = os.path.join(config["outdir"],"local_trees","{tree}.txt")
    shell:
        """
        clusterfunk focus -i {input.tree:q} \
        -o {output.collapsed_tree:q} \
        --metadata {input.metadata:q} \
        --index-column closest \
        --in-format nexus \
        --out-format newick \
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
