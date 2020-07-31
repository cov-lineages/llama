#!/usr/bin/env python3
import os
from matplotlib import font_manager as fm, rcParams
import baltic as bt

import re
import copy

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.patheffects as path_effects
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

import numpy as np
from scipy.special import binom
import math

import itertools
import requests

from io import StringIO as sio
from io import BytesIO as csio
from Bio import Phylo
from collections import defaultdict

import datetime as dt
from collections import Counter
from collections import defaultdict


thisdir = os.path.abspath(os.path.dirname(__file__))

def find_tallest_tree(input_dir):
    tree_heights = []
    
    for r,d,f in os.walk(input_dir):
        for fn in f:
            if fn.endswith(".tree"):
               
                tree_file = os.path.join(r, fn)
                tree = bt.loadNewick(tree_file,absoluteTime=False)
                tips = []
                
                for k in tree.Objects:
                    if k.branchType == 'leaf':
                        tips.append(k.name)
                
                tree_heights.append(tree.treeHeight)
    
    max_height = sorted(tree_heights, reverse=True)[0]
    return max_height

def display_name(tree, tree_name, tree_dir, full_taxon_dict, query_dict, label_fields):
    for k in tree.Objects:
        if k.branchType == 'leaf':
            name = k.name
            
            if "inserted" in name:
                collapsed_node_info = summarise_collapsed_node_for_label(tree_dir, name, tree_name, full_taxon_dict)
                k.traits["display"] = collapsed_node_info
            else:
                if name in full_taxon_dict:
                    taxon_obj = full_taxon_dict[name]
                
                    date = taxon_obj.sample_date
                    global_lineage = taxon_obj.global_lin
                    
                    k.traits["display"] = f"{name}|{date}|{global_lineage}"    

                    

                    if name in query_dict.keys():
                        if len(label_fields) > 0: 
                            for label_element in label_fields:
                                k.traits["display"] = k.traits["display"] + "|" + taxon_obj.attribute_dict[label_element]            
                
                else:
                    k.traits["display"] = name + "|" + "not in dict"


def find_colour_dict(query_dict, trait): 

    attribute_options = set()

    cmap = cm.get_cmap("viridis")

    if trait == "adm1": 
        colour_dict = {"Wales":"darkseagreen",
                "England":"indianred",
                "Scotland":"steelblue",
                "Northern_Ireland":"skyblue",
                "NA": "goldenrod"}
        return colour_dict

    else:
        for query in query_dict.values():
            attribute_options.add(query.attribute_dict[trait])
            
    if len(attribute_options) == 2:
        colour_dict = {list(attribute_options)[0]: "goldenrod",
                        list(attribute_options)[1]:"midnightblue"}
        return colour_dict

    else:
        #get the right number of colours, then loop through the set
        colour_dict = {}
        count = 0
        colors = cmap(np.linspace(0, 1, len(attribute_options)))
        for option in attribute_options:
            colour_dict[option] = colors[count]
            count += 1
    
        return colour_dict


def make_scaled_tree_without_legend(My_Tree, tree_name, tree_dir, num_tips, colour_dict_dict, colour_fields, label_fields, tallest_height,lineage, taxon_dict, query_dict):
#make colour_dict_dict optional argument
    display_name(My_Tree, tree_name, tree_dir, taxon_dict, query_dict, label_fields) 
    My_Tree.uncollapseSubtree()


    if num_tips < 10:
        #page_height = num_tips/2
        page_height = num_tips
    else:
        #page_height = num_tips/4 
        page_height = num_tips/2  

    offset = tallest_height - My_Tree.treeHeight
    space_offset = tallest_height/10
    absolute_x_axis_size = tallest_height+space_offset+space_offset + tallest_height #changed from /3 
    
    tipsize = 40
    c_func=lambda k: 'dimgrey' ## colour of branches
    l_func=lambda k: 'lightgrey' ## colour of branches
    s_func = lambda k: tipsize*5 if k.name in k.name in query_dict.keys() else tipsize
    z_func=lambda k: 100
    b_func=lambda k: 2.0 #branch width
    so_func=lambda k: tipsize*5 if k.name in k.name in query_dict.keys() else 0
    zo_func=lambda k: 99
    zb_func=lambda k: 98
    zt_func=lambda k: 97
    font_size_func = lambda k: 25 if k.name in k.name in query_dict.keys() else 15
    kwargs={'ha':'left','va':'center','size':12}

    if colour_fields != []:
        trait = colour_fields[0] #so always have the first trait as the first colour dot
        colour_dict = colour_dict_dict[trait]

        cn_func = lambda k: colour_dict[query_dict[k.name].attribute_dict[trait]] if k.name in query_dict.keys() else 'dimgrey'
        co_func=lambda k: colour_dict[query_dict[k.name].attribute_dict[trait]] if k.name in query_dict.keys() else 'dimgrey' 
        outline_colour_func = lambda k: colour_dict[query_dict[k.name].attribute_dict[trait]] if k.name in query_dict.keys() else 'dimgrey' 

    else:

        cn_func = lambda k: "goldenrod" if k.name in query_dict.keys() else 'dimgrey'
        co_func=lambda k: "goldenrod" if k.name in query_dict.keys() else 'dimgrey' 
        outline_colour_func = lambda k: "goldenrod" if k.name in query_dict.keys() else 'dimgrey' 

    x_attr=lambda k: k.height + offset
    y_attr=lambda k: k.y

    y_values = []
    for k in My_Tree.Objects:
        y_values.append(y_attr(k))
    min_y_prep = min(y_values)
    max_y_prep = max(y_values)
    vertical_spacer = 0.5 
    full_page = page_height + vertical_spacer + vertical_spacer
    min_y,max_y = min_y_prep-vertical_spacer,max_y_prep+vertical_spacer

    x_values = []
    for k in My_Tree.Objects:
        x_values.append(x_attr(k))
    max_x = max(x_values)
    
    
    fig2,ax2 = plt.subplots(figsize=(20,page_height),facecolor='w',frameon=False, dpi=100)
    

    My_Tree.plotTree(ax2, colour_function=c_func, x_attr=x_attr, y_attr=y_attr, branchWidth=b_func)
    My_Tree.plotPoints(ax2, x_attr=x_attr, colour_function=cn_func,y_attr=y_attr, size_function=s_func, outline_colour=outline_colour_func)
    My_Tree.plotPoints(ax2, x_attr=x_attr, colour_function=co_func, y_attr=y_attr, size_function=so_func, outline_colour=outline_colour_func)

    blob_dict = {}

    for k in My_Tree.Objects:
        
        if "display" in k.traits:
            name=k.traits["display"]

            x=x_attr(k)
            y=y_attr(k)
        
            height = My_Tree.treeHeight+offset
            text_start = tallest_height+space_offset+space_offset

            

            if len(colour_fields) > 1:
                
                division = (text_start - tallest_height)/(len(colour_fields))
                tip_point = tallest_height+space_offset

                if k.name in query_dict.keys():
                    
                    count = 0
                    
                    for trait in colour_fields[1:]:
                        
                        x_value = tip_point + count
                        count += division

                        option = query_dict[k.name].attribute_dict[trait]
                        colour_dict = colour_dict_dict[trait]
                        trait_blob = ax2.scatter(x_value, y, tipsize*5, color=colour_dict[option])  
                        
                        blob_dict[trait] = x_value

                    ax2.text(text_start+division, y, name, size=font_size_func(k), ha="left", va="center", fontweight="light")
                    if x != max_x:
                        ax2.plot([x+space_offset,tallest_height],[y,y],ls='--',lw=1,color=l_func(k))

                else:

                    ax2.text(text_start+division, y, name, size=font_size_func(k), ha="left", va="center", fontweight="light")
                    if x != max_x:
                        ax2.plot([x+space_offset,tallest_height],[y,y],ls='--',lw=1,color=l_func(k))


                for blob_x in blob_dict.values():

                    line_x = blob_x - (division/2)

                    ax2.plot([line_x,line_x],[min_y,max_y],ls='--',lw=3,color=l_func(k))
            
            
            else:
                ax2.text(text_start, y, name, size=font_size_func(k), ha="left", va="center", fontweight="ultralight")
                ax2.plot([x+space_offset,tallest_height+space_offset],[y,y],ls='--',lw=1,color=l_func(k))

    if len(colour_fields) > 1:

        blob_dict[colour_fields[0]] = tallest_height
        
        for trait, blob_x in blob_dict.items():

            y = max_y
            x = blob_x

            ax2.text(x,y,trait, rotation=45, size=15)


    ax2.spines['top'].set_visible(False) ## make axes invisible
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    
    ax2.set_xlim(-space_offset,absolute_x_axis_size)
    ax2.set_ylim(min_y,max_y)

    plt.yticks([])
    plt.xticks([])

    fig2.tight_layout()


def sort_trees_index(tree_dir):
    b_list = []
    d_list = []
    for r,d,f in os.walk(tree_dir):
        for thing in f:
            if thing.endswith("tree"):
                a = thing.split(".")[0]
                b = a.split("_")[1]
                b_list.append(int(b))
        
    c = sorted(b_list, key=int)
        
    return c

def make_all_of_the_trees(input_dir, taxon_dict, query_dict, colour_fields, label_fields, min_uk_taxa=3):

    tallest_height = find_tallest_tree(input_dir)

    too_tall_trees = []
    colour_dict_dict = defaultdict(dict)

    overall_df_dict = defaultdict(dict)

    overall_tree_count = 0
    
    lst = sort_trees_index(input_dir)

    for trait in colour_fields:
        colour_dict = find_colour_dict(query_dict, trait)
        colour_dict_dict[trait] = colour_dict

    for tree_number in lst:
        treename = "tree_" + str(tree_number)
        treefile = "local_" + str(tree_number) + ".tree"
        nodefile = "local_" + str(tree_number)
        num_taxa = 0

        tree = bt.loadNewick(input_dir + "/" + treefile, absoluteTime=False)

        old_node = tree.root
        new_node = bt.node()
        new_node.children.append(old_node)
        old_node.parent = new_node
        old_node.length=2.0
        new_node.height = 0
        new_node.y = old_node.y
        tree.root = new_node

        tree.Objects.append(new_node)

        tips = []
        
        for k in tree.Objects:
            if k.branchType == 'leaf':
                tips.append(k.name)
        
        if len(tips) < 1000:

            df_dict = summarise_node_table(input_dir, nodefile, taxon_dict)

            overall_df_dict[treename] = df_dict

            overall_tree_count += 1     
        
            make_scaled_tree_without_legend(tree, nodefile, input_dir, len(tips), colour_dict_dict, colour_fields, label_fields,tallest_height, tree_number, taxon_dict, query_dict)   
  
        else:
            too_tall_trees.append(tree_number)
            continue

    return too_tall_trees, overall_tree_count, overall_df_dict, colour_dict_dict

def summarise_collapsed_node_for_label(tree_dir, focal_node, focal_tree, full_tax_dict): 
    
    focal_tree_file = focal_tree + ".txt"

    with open(tree_dir + "/" + focal_tree_file) as f:
        next(f)
        for l in f:
            toks = l.strip("\n").split("\t")
            node_name = toks[0]
            members = toks[1]
        
            if node_name == focal_node:
                lineages = []
                
                member_list = members.split(",")
                number_nodes = str(len(member_list)) + " nodes"

                for tax in member_list:
                    if tax in full_tax_dict.keys():
                        taxon_obj = full_tax_dict[tax]
                        
                        lineages.append(taxon_obj.global_lin)
                    
                    else: #should always be in the full metadata now
                        print("tax missing from full metadata")
                    
                lineage_counts = Counter(lineages)

                most_common_lineages = []

                if len(lineage_counts) > 5:
                    
                    remaining = len(lineage_counts) - 5
                    
                    most_common_tups = lineage_counts.most_common(5)
                    for i in most_common_tups:
                        most_common_lineages.append(i[0])

                    pretty_lineages_prep = str(most_common_lineages).lstrip("[").rstrip("]").replace("'", "")
                    
                    if remaining == 1:
                        pretty_lineages = pretty_lineages_prep + " and " + str(remaining) + " other"
                    else:
                        pretty_lineages = pretty_lineages_prep + " and " + str(remaining) + " others"
                
                else:
                    pretty_lineages = str(list(lineage_counts.keys())).lstrip("[").rstrip("]").replace("'", "")


                node_number = node_name.lstrip("inserted_node")
                pretty_node_name = "Collapsed node " + node_number

                info = pretty_node_name + ": " + number_nodes + " in " + pretty_lineages

    return info

def summarise_node_table(tree_dir, focal_tree, full_tax_dict):

    focal_tree_file = focal_tree + ".txt"

    df_dict = defaultdict(list)

    with open(tree_dir + "/" + focal_tree_file) as f:
        next(f)
        for l in f:
            toks = l.strip("\n").split("\t")
            node_name = toks[0]
            members = toks[1]
        
            dates = []
            countries = []

            node_number = node_name.lstrip("inserted_node")
            
            member_list = members.split(",")

            for tax in member_list:
                if tax in full_tax_dict.keys():
                    taxon_obj = full_tax_dict[tax]
                
                    if taxon_obj.sample_date != "NA":
                        date_string = taxon_obj.sample_date
                        date = dt.datetime.strptime(date_string, "%Y-%m-%d").date()
                        dates.append(date)
                    
                    countries.append(taxon_obj.country)


            country_counts = Counter(countries)

            most_commons = country_counts.most_common(5)

            country_str = ""

            elem_count = 0

            for country, count in most_commons:
                elem_count += 1
                if elem_count == len(most_commons):
                    elem = country + " (" + str(count) + ")"
                    country_str += elem
                else:
                    elem = country + " (" + str(count) + "), "
                    country_str += elem
                

            min_date = str(min(dates))
            max_date = str(max(dates))

            size = len(member_list)

            df_dict["Node number"].append(node_number)
            df_dict["Number of sequences"].append(size)
            df_dict["Date range"].append(min_date + " to " + max_date)
            df_dict["Countries"].append(country_str)

    return df_dict

def make_legend(colour_dict):
    
    fig,ax = plt.subplots(figsize=(len(colour_dict)+1,1))

    plt.gca().set_aspect('equal', adjustable='box')
    plt.text
    
    x = 0
    for option in colour_dict.keys():
        circle = plt.Circle((x, 0.5), 0.05, color=colour_dict[option]) #((xloc, yloc), radius) relative to overall plot size
        ax.add_artist(circle)
        plt.text(x-0.1,0.3,option, fontsize=5)
        x += 1
        
        
    length = len(colour_dict)

    plt.xlim(-1,length)
    plt.ylim(0,1)

    ax.spines['top'].set_visible(False) ## make axes invisible
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)


    plt.yticks([])
    plt.xticks([])
    plt.show()

def describe_lineages(full_tax_dict):

    lineages_prep = defaultdict(list)
    lineages_present = defaultdict(dict)

    for tax in full_tax_dict.values():
        if tax.tree != "NA":
            key = tax.tree 
            lineages_prep[key].append(tax.global_lin)

    for tree, lineages in lineages_prep.items():
        counts = Counter(lineages)
        lineages_present[tree] = counts

    fig_count = 1
    tree_to_lin_fig = {}
    for tree, counts in lineages_present.items():
        
        fig, ax = plt.subplots(1,1, figsize=(5,5), dpi=100)

        x = list(counts.keys())
        y = list(counts.values())

        ax.bar(x,y)

        ax.set_xticklabels(x, rotation=45)
        ax.set_ylabel("Number of sequences")
        ax.set_xlabel("Global lineage")
        
        tree_to_lin_fig[tree] = fig_count
        fig_count += 1
    
    return tree_to_lin_fig
    


