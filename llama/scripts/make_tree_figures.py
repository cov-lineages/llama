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

def display_name(tree, tree_name, tree_dir, full_taxon_dict):
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
                
                else:
                    k.traits["display"] = name + "|" + "not in dict"


def make_scaled_tree_without_legend(My_Tree, tree_name, tree_dir, num_tips, tallest_height,lineage, taxon_dict, query_dict):

    display_name(My_Tree, tree_name, tree_dir, taxon_dict) 
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

    cn_func = lambda k: "goldenrod" if k.name in k.name in query_dict.keys() else 'dimgrey'
    co_func=lambda k: "goldenrod" if k.name in k.name in query_dict.keys() else 'dimgrey' 
    outline_colour_func = lambda k:"goldenrod" if k.name in k.name in query_dict.keys() else 'dimgrey' 

    x_attr=lambda k: k.height + offset
    y_attr=lambda k: k.y

    y_values = []
    for k in My_Tree.Objects:
        y_values.append(y_attr(k))
    
    y_values = sorted(y_values)
    vertical_spacer = 0.5 
    full_page = page_height + vertical_spacer + vertical_spacer
    min_y,max_y = y_values[0]-vertical_spacer,y_values[-1]+vertical_spacer
    
    
    fig2,ax2 = plt.subplots(figsize=(20,page_height),facecolor='w',frameon=False, dpi=100)
    

    My_Tree.plotTree(ax2, colour_function=c_func, x_attr=x_attr, y_attr=y_attr, branchWidth=b_func)
    My_Tree.plotPoints(ax2, x_attr=x_attr, colour_function=cn_func,y_attr=y_attr, size_function=s_func, outline_colour=outline_colour_func)
    My_Tree.plotPoints(ax2, x_attr=x_attr, colour_function=co_func, y_attr=y_attr, size_function=so_func, outline_colour=outline_colour_func)


    for k in My_Tree.Objects:
        
        if "display" in k.traits:
            name=k.traits["display"]

            x=x_attr(k)
            y=y_attr(k)
        
            height = My_Tree.treeHeight+offset
            
            ax2.text(tallest_height+space_offset+space_offset, y, name, size=font_size_func(k), ha="left", va="center", fontweight="ultralight")
            ax2.plot([x+space_offset,tallest_height+space_offset],[y,y],ls='--',lw=0.5,color=l_func(k))


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

def make_all_of_the_trees(input_dir, taxon_dict, query_dict, min_uk_taxa=3):

    tallest_height = find_tallest_tree(input_dir)

    too_tall_trees = []

    overall_df_dict = defaultdict(dict)

    overall_tree_count = 0
    
    lst = sort_trees_index(input_dir)

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
        
            make_scaled_tree_without_legend(tree, nodefile, input_dir, len(tips), tallest_height, tree_number, taxon_dict, query_dict)     
        else:
            too_tall_trees.append(tree_number)
            continue

    return too_tall_trees, overall_tree_count, overall_df_dict

def summarise_collapsed_node_for_label(tree_dir, focal_node, focal_tree, full_tax_dict): #could have lineages instead?
    
    focal_tree_file = focal_tree + ".txt"

    with open(tree_dir + "/" + focal_tree_file) as f:
        next(f)
        for l in f:
            toks = l.strip("\n").split("\t")
            node_name = toks[0]
            members = toks[1]
        
            if node_name == focal_node:
                countries = []
                
                member_list = members.split(",")
                number_nodes = str(len(member_list)) + " nodes"

                for tax in member_list:
                    if tax in full_tax_dict.keys():
                        taxon_obj = full_tax_dict[tax]
                        
                        countries.append(taxon_obj.country)
                    
                    else: #should always be in the full metadata now
                        print("tax missing from full metadata")
                    
                country_counts = Counter(countries)

                most_common_countries = []

                if len(country_counts) > 5:
                    
                    remaining = len(country_counts) - 5
                    
                    most_common_tups = country_counts.most_common(5)
                    for i in most_common_tups:
                        most_common_countries.append(i[0])

                    pretty_countries_prep = str(most_common_countries).lstrip("[").rstrip("]").replace("'", "")
                    
                    if remaining == 1:
                        pretty_countries = pretty_countries_prep + " and " + str(remaining) + " other"
                    else:
                        pretty_countries = pretty_countries_prep + " and " + str(remaining) + " others"
                
                else:
                    pretty_countries = str(list(country_counts.keys())).lstrip("[").rstrip("]").replace("'", "")


                node_number = node_name.lstrip("inserted_node")
                pretty_node_name = "Collapsed node " + node_number

                info = pretty_node_name + ": " + number_nodes + " in " + pretty_countries

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


def describe_tree_background(full_tax_dict, tree_dir, node_table):

    tree_lst = sort_trees_index(tree_dir)
    bar_to_tree = {}

    hidden_countries = defaultdict(list)

    figure_count = 0

    for fn in tree_lst:
        focal_tree = "local_" + str(fn)
        focal_tree_file = tree_dir + "/" + focal_tree + ".txt"
        pretty_focal = "Tree " + str(fn)

        collapsed_dict = defaultdict(list)

        with open(focal_tree_file) as f:
            next(f)
            for l in f:
                toks = l.strip("\n").split("\t")
                seqs = toks[1].split(",")

                node_number = toks[0].lstrip("inserted_node")
                new_name = "Collapsed node" + node_number
                collapsed_dict[new_name] = seqs
    
            ndes_country_counts = defaultdict(dict)
            nodes = []

            for nde, seqs in collapsed_dict.items():
                countries = []
                for i in seqs:
                    obj = full_tax_dict[i]
                    countries.append(obj.country)

                    country_counts = Counter(countries)
                                
                if len(country_counts) > 10:
                    keep_countries = dict(country_counts.most_common(10))

                    hidden_countries[focal_tree].append(nde)
                    
                else:
                    keep_countries = country_counts
                    
                if len(country_counts) > 1:
                    
                    ndes_country_counts[nde] = keep_countries
                    nodes.append(nde)
                            
            if len(ndes_country_counts) > 1:
                
                figure_count += 1

                plt.rc('ytick', labelsize=5)
                
                count = 0

                rows = math.ceil(len(ndes_country_counts)/5)
                

                if rows == 1:
                    fig, axs = plt.subplots(rows,5, figsize=(10,2)) 

                    fig.tight_layout()
                    count = 0      
                    for nde, country_counts in ndes_country_counts.items():

                        x = country_counts.keys()
                        y = country_counts.values()

                        # print("1 counts" + str(count))
                        axs[count].bar(x,y, color="goldenrod")
                        axs[count].set_title(nde, size=8)
                        axs[count].set_xticklabels(x,rotation=90, size=5)
                        #axs[count].set_yticklabels(size=5)
                        
                        count += 1

                 
                    fig.suptitle(pretty_focal,y=1.1,x=0.05, size=10)
                
                else:
                    fig, axs = plt.subplots(rows,5,figsize=(10,10))
                    fig.subplots_adjust(hspace=1.0, wspace=0.7)
                    # fig.tight_layout()
                    
                    
                    for nrow in range(0,rows):
                        for i in range(0,5):
                            try:
                                relevant_nde = nodes[(nrow*5) + i]
                                
                                x = ndes_country_counts[relevant_nde].keys()
                                y = ndes_country_counts[relevant_nde].values()

                                axs[nrow][i].bar(x,y, color="goldenrod")
                                axs[nrow][i].set_title(relevant_nde, size=8)
                                axs[nrow][i].set_xticklabels(x,rotation=70, size=5)
                                # axs[nrow][i].set_yticklabels(y, size=5)
                            except IndexError:
                                continue

               
                    fig.suptitle(pretty_focal,y=0.95,x=0.1, size=10)

                if len(ndes_country_counts) != rows*5:
                    number_empty_ones = rows*5 - len(ndes_country_counts)
                    for_removal = [i for i in range((rows*5-number_empty_ones),rows*5)]

                    for j in for_removal:
                         fig.delaxes(axs.flatten()[j])

                


                    
            elif len(ndes_country_counts) == 1:
                
                figure_count += 1
                plt.figure(figsize=(2,2))

                for nde, country_counts in ndes_country_counts.items():
                    
                    x = country_counts.keys()
                    y = country_counts.values()

                    plt.bar(x,y, color="goldenrod")
                    # plt.title(nde)
                    plt.xticks(size=5, rotation=90)
                    plt.yticks(size=5)

                    plt.title(pretty_focal + ": " + nde, size=5)


            bar_to_tree[figure_count] = focal_tree


    return figure_count, bar_to_tree

