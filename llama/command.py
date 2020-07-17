#!/usr/bin/env python3
from llama import __version__
import setuptools
import argparse
import os.path
import snakemake
import sys
from tempfile import gettempdir
import tempfile
import pprint
import json
import csv
import os
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
    usage='''llama <query> [options]''')

    parser.add_argument('query',help="Input csv file with minimally `name` as a column header. Alternatively, `--index-column` can specifiy a column name other than `name`")
    parser.add_argument('-i',"--id-string", action="store_true",help="Indicates the input is a comma-separated id string with one or more query ids. Example: `EDB3588,EDB3589`.", dest="ids")
    parser.add_argument('--fasta', action="store",help="Optional fasta query.", dest="fasta")
    parser.add_argument('-o','--outdir', action="store",help="Output directory. Default: current working directory")
    parser.add_argument('--datadir', action="store",help="Local directory that contains the data files")
    parser.add_argument('--index-column', action="store",help="Input csv column to match in database. Default: name", dest="index_column",default="name")
    parser.add_argument('--search-field', action="store",help="Column in database to match with input csv. Default: covv_accession_id", dest="search_field",default="covv_accession_id")
    parser.add_argument('--distance', action="store",help="Extraction from large tree radius. Default: 2", dest="distance",default=2)
    parser.add_argument('-n', '--dry-run', action='store_true',help="Go through the motions but don't actually run")
    parser.add_argument('--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: $TMPDIR")
    parser.add_argument("--no-temp",action="store_true",help="Output all intermediate files, for dev purposes.")
    parser.add_argument('--collapse-threshold', action='store',type=int,help="Minimum number of nodes to collapse on. Default: 1", dest="threshold", default=1)
    parser.add_argument('-t', '--threads', action='store',type=int,help="Number of threads")
    parser.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    parser.add_argument('--max-ambig', action="store", default=0.5, type=float,help="Maximum proportion of Ns allowed to attempt analysis. Default: 0.5",dest="maxambig")
    parser.add_argument('--min-length', action="store", default=10000, type=int,help="Minimum query length allowed to attempt analysis. Default: 10000",dest="minlen")
    parser.add_argument("-v","--version", action='version', version=f"llama {__version__}")

    # Exit with help menu if no args supplied
    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    # find the master Snakefile
    snakefile = os.path.join(thisdir, 'scripts','Snakefile')
    if not os.path.exists(snakefile):
        sys.stderr.write('Error: cannot find Snakefile at {}\n Check installation'.format(snakefile))
        sys.exit(-1)
    
    # find the query fasta
    if args.fasta:
        fasta = os.path.join(cwd, args.fasta)
        if not os.path.exists(fasta):
            sys.stderr.write('Error: cannot find fasta query at {}\n'.format(fasta))
            sys.exit(-1)
        else:
            print(f"Input fasta file: {fasta}")
    else:
        fasta = ""

    # default output dir
    outdir = ''
    if args.outdir:
        rel_outdir = args.outdir #for report weaving
        outdir = os.path.join(cwd, args.outdir)
        
        if not os.path.exists(outdir):
            os.mkdir(outdir)
    else:
        timestamp = str(datetime.now().isoformat(timespec='milliseconds')).replace(":","").replace(".","").replace("T","-")
        outdir = os.path.join(cwd, timestamp)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        rel_outdir = os.path.join(".",timestamp)
    
    print(f"Output files will be written to {outdir}\n")

    # specifying temp directory
    tempdir = ''
    if args.tempdir:
        to_be_dir = os.path.join(cwd, args.tempdir)
        if not os.path.exists(to_be_dir):
            os.mkdir(to_be_dir)
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=to_be_dir)
        tempdir = temporary_directory.name
    else:
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=None)
        tempdir = temporary_directory.name

    # if no temp, just write everything to outdir
    if args.no_temp:
        print(f"--no-temp: All intermediate files will be written to {outdir}")
        tempdir = outdir

    # find the query csv, or string of ids
    query = os.path.join(cwd, args.query)
    if not os.path.exists(query):
        if args.ids:
            id_list = args.query.split(",")
            query = os.path.join(tempdir, "query.csv")
            with open(query,"w") as fw:
                fw.write(f"{args.index_column}\n")
                for i in id_list:
                    fw.write(i+'\n')
        else:
            sys.stderr.write(f"Error: cannot find query file at {query}\n Check if the file exists, or if you're inputting an id string (e.g. EPI12345,EPI23456), please use in conjunction with the `--id-string` flag\n.")
            sys.exit(-1)
    else:
        print(f"Input file: {query}")

    # parse the input csv, check col headers
    queries = []
    with open(query, newline="") as f:
        reader = csv.DictReader(f)
        column_names = reader.fieldnames
        if args.index_column not in column_names:
            sys.stderr.write(f"Error: Input file missing header field {args.index_column}\nEither specifiy `--index-column` or supply a column with header `name`\n")
            sys.exit(-1)
                    
        print("Input querys to process:")
        for row in reader:
            queries.append(row[args.index_column])
            
            print(row[args.index_column])
        print(f"Total: {len(queries)}")
    print('\n')
    # how many threads to pass
    if args.threads:
        threads = args.threads
    else:
        threads = 1
    print(f"Number of threads: {threads}\n")

    # create the config dict to pass through to the snakemake file
    config = {
        "query":query,
        "outdir":outdir,
        "tempdir":tempdir,
        "trim_start":265,   # where to pad to using datafunk
        "trim_end":29674,   # where to pad after using datafunk
        "fasta":fasta,
        "rel_outdir":rel_outdir,
        "search_field":args.search_field,
        "index_column":args.index_column,
        "force":"True"
        }

    # find the data files
    data_dir = ""
    if args.datadir:
        data_dir = os.path.join(cwd, args.datadir)
        metadata,seqs,tree = ("","","")
        
        seqs = os.path.join(data_dir,"alignment.fasta")
        
        metadata = os.path.join(data_dir,"metadata.csv")

        tree = os.path.join(data_dir,"global.tree")

        if not os.path.isfile(seqs) or not os.path.isfile(metadata) or not os.path.isfile(tree):
            sys.stderr.write(f"""Error: cannot find correct data files at {data_dir}\nThe directory should contain the following files:\n\
    - alignment.fasta\n\
    - global.tree\n\
    - metadata.csv\n""")
            sys.exit(-1)
        else:
            config["seqs"] = seqs
            config["metadata"] = metadata
            config["tree"] = tree

            print("Found data:")
            print("    -",seqs)
            print("    -",metadata)
            print("    -",tree,"\n")
    else:
        print("No data directory specified, please specify where to find the data files\n")
        sys.exit(-1)

    """ 
    QC steps:
    1) check csv header
    2) check fasta file N content
    3) write a file that contains just the seqs to run
    """
    # run qc on the input sequence file
    if args.fasta:
        do_not_run = []
        run = []
        for record in SeqIO.parse(args.fasta, "fasta"):
            if len(record) <args.minlen:
                record.description = record.description + f" fail=seq_len:{len(record)}"
                do_not_run.append(record)
                print(f"    - {record.id}\tsequence too short: Sequence length {len(record)}")
            else:
                num_N = str(record.seq).upper().count("N")
                prop_N = round((num_N)/len(record.seq), 2)
                if prop_N > args.maxambig: 
                    record.description = record.description + f" fail=N_content:{prop_N}"
                    do_not_run.append(record)
                    print(f"    - {record.id}\thas an N content of {prop_N}")
                else:
                    run.append(record)

        post_qc_query = os.path.join(outdir, 'query.post_qc.fasta')
        with open(post_qc_query,"w") as fw:
            SeqIO.write(run, fw, "fasta")
        qc_fail = os.path.join(outdir,'query.failed_qc.csv')
        with open(qc_fail,"w") as fw:
            fw.write("name,reason_for_failure\n")
            for record in do_not_run:
                desc = record.description.split(" ")
                for i in desc:
                    if i.startswith("fail="):
                        fw.write(f"{record.id},{i}\n")

        config["post_qc_query"] = post_qc_query
        config["qc_fail"] = qc_fail
    else:
        config["post_qc_query"] = ""
        config["qc_fail"] = ""


    # accessing package data and adding to config dict
    reference_fasta = pkg_resources.resource_filename('llama', 'data/reference.fasta')
    config["reference_fasta"] = reference_fasta

    if args.distance:
        try:
            distance = int(args.distance) 
            config["distance"] = args.distance
        except:
            sys.stderr.write('Error: distance must be an integer\n')
            sys.exit(-1)
    else:
        config["distance"] = "1"

    
    if args.threshold:
        try:
            threshold = int(args.threshold)
            config["threshold"] = args.threshold
        except:
            sys.stderr.write('Error: threshold must be an integer\n')
            sys.exit(-1)
    else:
        config["threshold"] = "1"

    # don't run in quiet mode if verbose specified
    if args.verbose:
        quiet_mode = False
        config["quiet_mode"]="False"
    else:
        quiet_mode = True
        config["quiet_mode"]="True"

    status = snakemake.snakemake(snakefile, printshellcmds=True,
                                 dryrun=args.dry_run, forceall=True,force_incomplete=True,workdir=tempdir,
                                 config=config, cores=threads,lock=False,quiet=quiet_mode
                                 )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1

if __name__ == '__main__':
    main()