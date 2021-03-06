#!/usr/bin/env python3
import os
from pweave import weave
import argparse
import shutil
import sys


thisdir = os.path.abspath(os.path.dirname(__file__))

def make_report(metadata, input_csv, filtered_metadata, outfile, outdir, full_outdir, treedir, figdir, report_template, failed_seqs, no_seq, input_column, data_column, label_fields, colour_fields, node_summary):

    name_stem = ".".join(outfile.split(".")[:-1])

    tree_name_stems = []
    for r,d,f in os.walk(treedir):
        for fn in f:
            if fn.endswith(".tree"):
                tree_name_stems.append(fn.split("_")[0])
    tree_name_stems = list(set(tree_name_stems))
    if len(tree_name_stems) > 1:
        sys.stderr.write("Error: Multiple tree names found")
        sys.exit(-1)
    elif len(tree_name_stems) == 0:
        sys.stderr.write("Error: No trees found in tree directory. Note, tree name much end with .tree")
        sys.exit(-1)
    else:
        tree_name_stem = tree_name_stems[0]

    with open(outfile, 'w') as pmd_file:
    
        md_template = report_template
        summary_dir = os.path.join(outdir, "summary_files")
 
        change_line_dict = {
                            "output_directory": f'output_directory = "{outdir}"\n',
                            "full_outdir": f'full_outdir = "{full_outdir}"\n',
                            "name_stem_input": f'name_stem_input = "{name_stem}"\n',
                            "full_metadata_file": f'full_metadata_file = "{metadata}"\n',
                            "filtered_db_metadata": f'filtered_db_metadata = "{filtered_metadata}"\n',
                            "input_csv": f'input_csv = "{input_csv}"\n',
                            "input_directory": f'input_directory = "{treedir}"\n',
                            "figdir": f'figdir = "{figdir}"\n',
                            "tree_dir": f'tree_dir = "{treedir}"\n',
                            "tree_name_stem": f'tree_name_stem = "{tree_name_stem}"\n',
                            "summary_dir": f'summary_dir = "{summary_dir}"\n',
                            "QC_fail_file": f'QC_fail_file = "{failed_seqs}"\n',
                            "missing_seq_file": f'missing_seq_file = "{no_seq}"\n',
                            "input_name_column": f'input_name_column = "{input_column}"\n',
                            "database_name_column": f'database_name_column = "{data_column}"\n',
                            "colour_fields": f'colour_fields_input = "{colour_fields}"\n',
                            "label_fields": f'label_fields_input = "{label_fields}"\n',
                            "node_summary_option": f'node_summary_option = "{node_summary}"\n'
                            }
        with open(md_template) as f:
            for l in f:
                if "##CHANGE" in l:
                    for key in change_line_dict:
                        if key in l:
                            new_l = change_line_dict[key]
                else:
                    new_l = l

                pmd_file.write(new_l)

    weave(outfile, doctype = "pandoc", figdir=figdir)

def main():
    parser = argparse.ArgumentParser(description="Report generator script")
    parser.add_argument("-i","--input-csv", required=False, help="path to input file",dest="input_csv")

    parser.add_argument("--filtered-metadata", required=False, help="path to combined metadata file",dest="filtered_metadata")
    parser.add_argument("--metadata", required=True, help="path to full metadata file",dest="metadata")
    
    parser.add_argument("--failed-seqs", required=False, default="", help="csv of seqs that fail qc and the reason why",dest="failed_seqs")
    parser.add_argument("--no-seq-provided", required=False, default="", help="file of seqs that weren't in cog and didn't have a sequence provided",dest="no_seq")
    
    parser.add_argument("-t","--treedir", required=False, default="", help="path to tree directory",dest="treedir")
    parser.add_argument("--report-template", help="report template file",dest="report_template")

    parser.add_argument("-o","--outfile", default="llama_report.pmd", help="output name stem as a string",dest="outfile")
    parser.add_argument("--outdir", help="output directory",dest="outdir")
    parser.add_argument("--full-outdir", help="full output directory",dest="full_outdir")

    parser.add_argument("--figdir", help="output directory",dest="figdir")

    parser.add_argument('--input-column', action="store",help="Column in input csv file to match with database. Default: name", dest="input_column",default="name")
    parser.add_argument('--data-column', action="store",help="Column in database to match with input csv file. Default: sequence_name", dest="data_column",default="sequence_name")

    parser.add_argument("--colour-fields", action="store", help="fields to colour tree tips by", dest="colour_fields")
    parser.add_argument("--label-fields", action="store", help="fields to label tree tips by", dest="label_fields")
    parser.add_argument("--node-summary", action="store", help="field to summarise collapsed nodes by. Default=lineage", dest="node_summary")


    args = parser.parse_args()

    make_report(args.metadata, args.input_csv, args.filtered_metadata, args.outfile, args.outdir, args.full_outdir, args.treedir, args.figdir, args.report_template, args.failed_seqs,args.no_seq, args.input_column, args.data_column,  args.label_fields, args.colour_fields, args.node_summary)


if __name__ == "__main__":
    main()