#!/usr/bin/env python3
from collections import defaultdict
import pandas as pd
import csv
from tabulate import tabulate
import baltic as bt
import os
import datetime as dt
import math


class taxon(): #might want country in here, not sure yet

    def __init__(self, name, global_lin):

        self.name = name

        self.sample_date = "NA"
        
        if global_lin == "":
            self.global_lin = "NA"
        else:
            self.global_lin = global_lin
       
        self.in_db = False

        self.tree = "NA"


    

def parse_filtered_metadata(metadata_file, tip_to_tree):
    
    query_dict = {}
    query_id_dict = {}
    present_lins = set()

    tree_to_tip = defaultdict(list)

    with open(metadata_file, "r") as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            glob_lin = sequence['lineage']

            query_id = sequence['query_id']
            query_name = sequence['query']
            closest_name = sequence["closest"]

            sample_date = sequence["sample_date"]

            new_taxon = taxon(query_name, glob_lin)

            new_taxon.query_id = query_id

            if query_name == closest_name: #if it's in the database, get its sample date
                new_taxon.in_db = True
                new_taxon.sample_date = sample_date
                new_taxon.closest = "NA"

            else:
                new_taxon.closest = closest_name

            present_lins.append(glob_lin)
            
            relevant_tree = tip_to_tree[query_name]
            new_taxon.tree = relevant_tree

            tree_to_tip[relevant_tree].append(new_taxon)
           
            query_dict[query_name] = new_taxon
            query_id_dict[query_id] = new_taxon
            
    return query_dict, query_id_dict, present_lins, tree_to_tip

def parse_input_csv(input_csv, query_id_dict, name_column):

    new_query_dict = {}

    with open(input_csv, 'r') as f:
        reader = csv.DictReader(f)
        col_name_prep = next(reader)
        col_names = list(col_name_prep.keys())

    with open(input_csv, 'r') as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            
            name = sequence[name_column] 

            if name in query_id_dict.keys():
                taxon = query_id_dict[name]

                if "sample_date" in col_names: #if it's not in database but date is provided (if it's in the database, it will already have been assigned a sample date.)
                    if sequence["sample_date"] != "":
                        taxon.sample_date = sequence["sample_date"]

                if "global_lineage" in col_names:
                    if taxon.global_lin == "NA" and sequence["global_lineage"] != "":
                        taxon.global_lin = sequence["global_lineage"]

                
                new_query_dict[taxon.name] = taxon
            
                
    return new_query_dict 

def parse_tree_tips(tree_dir):

    tips = []
    tip_to_tree = {}

    for fn in os.listdir(tree_dir):
        if fn.endswith("tree"):
            tree_name = fn.split(".")[0]
            tree = bt.loadNewick(tree_dir + "/" + fn, absoluteTime=False)
            for k in tree.Objects:
                if k.branchType == 'leaf' and "inserted" not in k.name:
                    tips.append(k.name)
                    tip_to_tree[k.name] = tree_name

        elif fn.endswith(".txt"):
            with open(tree_dir + "/" + fn) as f:
                for l in f:
                    tip_string = l.strip("\n").split("\t")[1]
                    tip_list = tip_string.split(",")
                    tips.extend(tip_list)

    return tips, tip_to_tree

def parse_full_metadata(query_dict, full_metadata, present_in_tree, database_name_column):

    full_tax_dict = query_dict.copy()

    with open(full_metadata, 'r') as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            seq_name = sequence[database_name_column] 

            date = sequence["sample_date"]

            country = sequence["country"]

            glob_lin = sequence["lineage"]

            if seq_name in present_in_tree and seq_name not in query_dict.keys():
                new_taxon = taxon(seq_name, glob_lin)
                if date == "":
                    date = "NA"
                
                new_taxon.sample_date = date

                full_tax_dict[seq_name] = new_taxon
                                    
    return full_tax_dict
    

def make_initial_table(query_dict):

    df_dict = defaultdict(list)

    for query in query_dict.values():
        
        df_dict["Query ID"].append(query.query_id.replace("|","\|"))
        
        if query.in_db: 
            df_dict["Sequence name in Tree"].append(query.name)
        else:
            df_dict["Sequence name in Tree"].append("NA")
        

        df_dict["Sample date"].append(query.sample_date)

        df_dict["Closest sequence in Tree"].append(query.closest)
        
        df_dict["Global lineage"].append(query.global_lin)
        
        if query.tree != "NA":
            tree_number = query.tree.split("_")[1]
            pretty_tree = "Tree " + str(tree_number)
            df_dict["Tree"].append(pretty_tree)
        else:
            df_dict["Tree"].append("NA") #this should never happen, it's more error catching

    df = pd.DataFrame(df_dict)

    df.set_index("Query ID", inplace=True)

    return df

def investigate_QC_fails(QC_file):

    fail_dict = {}

    with open(QC_file) as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            name = sequence["name"]
            reason = sequence["reason_for_failure"]

            if "seq_len" in reason:
                length = reason.split(":")[1]
                final_reason = "Sequence too short: only " + length + " bases."
            elif "N_content" in reason:
                n_content = reason.split(":")[1]
                final_reason = "Sequence has too many Ns: " + str(float(round(float(n_content)*100))) + "\% of bases"

            fail_dict[name] = final_reason


    return fail_dict

def print_missing_seqs(missing_seqs_file):
    
    failed_names = []

    with open(missing_seqs_file) as f:
        for l in f:
            name = l.strip("\n").split(",")[0]
            
            failed_names.append(name)

    return failed_names

