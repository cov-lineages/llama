#!/usr/bin/env python3
import os
from pweave import *
import argparse
import shutil


thisdir = os.path.abspath(os.path.dirname(__file__))

def make_report(metadata, input_csv, filtered_metadata, outfile, outdir, treedir, figdir, report_template, failed_seqs, no_seq):

    name_stem = ".".join(outfile.split(".")[:-1])
    with open(outfile, 'w') as pmd_file:
    
        md_template = report_template
        summary_dir = os.path.join(outdir, "summary_files")
 
        change_line_dict = {
                            "output_directory": f'output_directory = "{outdir}"\n',
                            "name_stem_input": f'name_stem_input = "{name_stem}"\n',
                            "full_metadata_file": f'full_metadata_file = "{metadata}"\n',
                            "filtered_metadata": f'filtered_metadata = "{filtered_metadata}"\n',
                            "input_csv": f'input_csv = "{input_csv}"\n',
                            "input_directory": f'input_directory = "{treedir}"\n',
                            "figdir": f'figdir = "{figdir}"\n',
                            "tree_dir": f'tree_dir = "{treedir}"\n',
                            "summary_dir": f'summary_dir = "{summary_dir}"\n',
                            "QC_fail_file": f'QC_fail_file = "{failed_seqs}"\n',
                            "missing_seq_file": f'missing_seq_file = "{no_seq}"\n'
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
    parser.add_argument("--figdir", help="output directory",dest="figdir")

    args = parser.parse_args()

    make_report(args.metadata, args.input_csv, args.filtered_metadata, args.outfile, args.outdir, args.treedir, args.figdir, args.report_template, args.failed_seqs,args.no_seq)


if __name__ == "__main__":
    main()