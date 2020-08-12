#!/usr/bin/env python3
from llama import __version__
import setuptools
import argparse
import os.path
import snakemake
import sys
from tempfile import gettempdir
import tempfile
import input_qc_functions as qcfunk
import pprint
import json
import csv
import os
import datetime 
from datetime import datetime
from Bio import SeqIO

import pkg_resources
from . import _program

"""
Need query_csv, metadata, fasta (opt)
"""


thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(prog = _program, 
    description='llama: Local Lineage And Monophyly Assessment', 
    usage='''llama -i <input.csv> -d <path/to/data> [options]''')

    parser.add_argument('-i','--input',help="Input csv file with minimally `name` as a column header. Alternatively, `--input-column` can specifiy a column name other than `name`",dest="query")
    parser.add_argument('-fm','--from-metadata',nargs='*', dest="from_metadata",help="Generate a query from the metadata file supplied. Define a search that will be used to pull out sequences of interest from the large phylogeny. E.g. -fm country=Ireland sample_date=2020-03-01:2020-04-01")
    parser.add_argument('-f','--fasta', action="store",help="Optional fasta query. Fasta sequence names must exactly match those in your input query.", dest="fasta")

    parser.add_argument('-a',"--align-sequences", action="store_true",help="Just align sequences.", dest="align")
    parser.add_argument('-s','--seqs', action="store",help="Sequence file containing sequences used to create the tree. For use in combination with the `--align-sequences` option.", dest="seqs")
    parser.add_argument('-ns','--no-seqs', action="store_true",help="Alignment not available. Note, to work, all queries must already be in global tree.", dest="no_seqs")
    
    parser.add_argument("-r","--report",action="store_true",help="Generate markdown report of input queries and their local trees")
    parser.add_argument('--colour-fields', action="store",help="Comma separated string of fields to colour by in the report.", dest="colour_fields")
    parser.add_argument('--label-fields', action="store", help="Comma separated string of fields to add to tree report labels.", dest="label_fields")

    parser.add_argument("--id-string", action="store_true",help="Indicates the input is a comma-separated id string with one or more query ids. Example: `EDB3588,EDB3589`.", dest="ids")
    parser.add_argument('-o','--outdir', action="store",help="Output directory. Default: current working directory")
    parser.add_argument('-d','--datadir', action="store",help="Local directory that contains the data files")
    parser.add_argument('--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: $TMPDIR")
    parser.add_argument("--no-temp",action="store_true",help="Output all intermediate files, for dev purposes.")
    
    parser.add_argument('--input-column', action="store",help="Column in input csv file to match with database. Default: name", dest="input_column",default="name")
    parser.add_argument('--data-column', action="store",help="Column in database to match with input csv file. Default: sequence_name", dest="data_column",default="sequence_name")
    
    parser.add_argument('--distance', action="store",type=int,help="Extraction from large tree radius. Default: 2", dest="distance",default=2)
    parser.add_argument('--collapse-threshold', action='store',type=int,help="Minimum number of nodes to collapse on. Default: 1", dest="threshold", default=1)
    parser.add_argument('--lineage-representatives', action="store_true",help="Include a selection of representative sequences from lineages present in the local tree. Default: False", dest="lineage_representatives")
    parser.add_argument('--number-of-representatives', action="store",type=int,help="How many representative sequeneces per lineage to keep in the collapsed tree. Default: 5", default=5, dest="number_of_representatives")

    parser.add_argument('--max-ambig', action="store", default=0.5, type=float,help="Maximum proportion of Ns allowed to attempt analysis. Default: 0.5",dest="maxambig")
    parser.add_argument('--min-length', action="store", default=10000, type=int,help="Minimum query length allowed to attempt analysis. Default: 10000",dest="minlen")
    
    parser.add_argument('-n', '--dry-run', action='store_true',help="Go through the motions but don't actually run")
    parser.add_argument('-t', '--threads', action='store',type=int,help="Number of threads")
    parser.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    parser.add_argument("--outgroup",action="store",help="Optional outgroup sequence to root local subtrees. Default an anonymised sequence that is at the base of the global SARS-CoV-2 phylogeny.")
    parser.add_argument("-v","--version", action='version', version=f"llama {__version__}")

    # Exit with help menu if no args supplied
    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    # find the Snakefile (either master, no-seq snakefile or alignment snakefile)
    snakefile = qcfunk.get_snakefile(args.no_seqs,args.align,thisdir)
    
    # find the query fasta
    fasta = qcfunk.get_query_fasta(args.fasta,args.no_seqs,cwd)

    # default output dir
    outdir,rel_outdir = qcfunk.get_outdir(args.outdir,cwd)

    print(f"Output files will be written to {outdir}\n")

    
    # if no temp, just write everything to outdir
    if args.no_temp:
        print(f"--no-temp: All intermediate files will be written to {outdir}")
        tempdir = outdir
    else:
        tempdir = qcfunk.get_temp_dir(args.tempdir, cwd)

    # how many threads to pass
    if args.threads:
        threads = args.threads
    else:
        threads = 1
    print(f"Number of threads: {threads}\n")

    # create the config dict to pass through to the snakemake file
    config = {
        "outdir":outdir,
        "tempdir":tempdir,
        "trim_start":265,   # where to pad to using datafunk
        "trim_end":29674,   # where to pad after using datafunk
        "fasta":fasta,
        "rel_outdir":rel_outdir,
        "data_column":args.data_column,
        "force":"True",
        "threshold":args.threshold,
        "number_of_representatives":args.number_of_representatives,
        "distance":args.distance
        }

    if args.lineage_representatives:
        config["lineage_representatives"]=args.lineage_representatives
    else:
        config["lineage_representatives"]=False

    # find the data files
    # data_dir = ""
    metadata,seqs,tree = ("","","")
    if not args.align:
        if args.datadir:
            metadata,seqs,tree = qcfunk.check_data_dir(args.datadir,args.no_seqs,cwd,config)
        else:
            sys.stderr.write(qcfunk.cyan("No data directory specified, please specify where to find the data files\n"))
            sys.exit(-1)

    elif args.align:
        seqs = qcfunk.get_seqs_for_aln(args.seqs,cwd)
        config["seqs"] = seqs

    if not args.align:
        # parse the input db, check col headers
        with open(metadata, newline="") as f:
            reader = csv.DictReader(f)
            column_names = reader.fieldnames
            if args.data_column not in column_names:
                sys.stderr.write(qcfunk.cyan(f"Error: Metadata file missing header field {args.data_column}\nEither specifiy `--data-column` or supply a column with header `sequence_name`\n"))
                sys.exit(-1)

        if args.from_metadata and args.query:
            sys.stderr.write(qcfunk.cyan(f"Error: Please provide either an input query file (`-i`) or define some search criteria from the metadata (`-fm`)\n"))
            sys.exit(-1)
        
        elif args.from_metadata:
            query = qcfunk.parse_from_metadata_arg(metadata, args.from_metadata, args.data_column, config)

        elif args.query:
            # find the query csv, or string of ids
            query = os.path.join(cwd, args.query)
            
            if not os.path.exists(query):
                if args.ids:
                    id_list = args.query.split(",")
                    query = os.path.join(tempdir, "query.csv")
                    with open(query,"w") as fw:
                        fw.write(f"{args.input_column}\n")
                        for i in id_list:
                            fw.write(i+'\n')
                else:
                    sys.stderr.write(qcfunk.cyan(f"Error: cannot find query file at {query}\nCheck if the file exists, or if you're inputting an id string (e.g. EPI12345,EPI23456), please use in conjunction with the `--id-string` flag\n."))
                    sys.exit(-1)
            else:
                print(f"Input file: {query}")
                
            config["query"] = query
            config["input_column"] = args.input_column
        else:
            sys.stderr.write(qcfunk.cyan(f"Error: please input a query (`-i`) or define a search (`-fm`)\n"))
            sys.exit(-1)
        # parse the input csv, check col headers

        if args.query:
            input_column = args.input_column
        else:
            input_column = args.data_column

        qcfunk.check_label_and_colour_fields(query, args.query, args.colour_fields, args.label_fields, input_column, config)


    if args.fasta:
        """ 
        QC steps:
        1) check csv header
        2) check fasta file N content
        3) write a file that contains just the seqs to run
        """
        qcfunk.input_file_qc(args.fasta,args.minlen,args.maxambig,config)
    else:
        config["post_qc_query"] = ""
        config["qc_fail"] = ""

    # accessing package data and adding to config dict
    qcfunk.get_outgroup_sequence(args.outgroup, cwd, config)

    if args.report:
        config["report"] = True
    else:
        config["report"] = False
    
    config["report_template"] =  os.path.join(thisdir, 'scripts','report_template.pmd')
    footer_fig = pkg_resources.resource_filename('llama', 'data/footer.png')
    config["footer"] = footer_fig

    # don't run in quiet mode if verbose specified
    if args.verbose:
        quiet_mode = False
        config["quiet_mode"]=False
    else:
        quiet_mode = True
        config["quiet_mode"]=True

    status = snakemake.snakemake(snakefile, printshellcmds=True,
                                 dryrun=args.dry_run, forceall=True,force_incomplete=True,workdir=tempdir,
                                 config=config, cores=threads,lock=False,quiet=quiet_mode
                                 )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1

if __name__ == '__main__':
    main()