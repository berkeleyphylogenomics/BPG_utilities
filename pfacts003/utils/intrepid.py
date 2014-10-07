#!/usr/bin/env python2.7

from Bio.PDB.Polypeptide import one_to_three

import os

def parse_intrepid_output(directory):
    """
    This function just parses the two intrepid outputs and returns the information from them necessary to populate
    the functional site prediction results page.
    TODO:
    Make this parser more general...
    """
    score_file = open(os.path.join(directory, "output.aux"), "r")
    rank_file = open(os.path.join(directory, "output.rank"), "r")
    rank_lines = rank_file.readlines()[1:]
    score_lines = score_file.readlines()[1:]
    table_row = []
    cons_js_score_array = []
    residue_name_array = []
    cons_js_rank_array = []
    for (index, line) in enumerate(rank_lines):
        (Position, ResSeq, Chain, Residue, ICode, cons_importance_map, cons_js_map, cons_re_map, global_map, global_js_map, global_re_map) = score_lines[index].replace("\n","").split("|")
        (cons_importance_map_rank, cons_js_map_rank, cons_re_map_rank, global_map_rank, global_js_map_rank, global_re_map_rank) = rank_lines[index].replace("\n","").split("|")[5:]
        residue_name_array.append(Residue + " (%s)" % one_to_three(Residue))
        cons_js_score_array.append(float(cons_js_map))
        table_row.append([int(Position), Residue + " (%s)" % one_to_three(Residue), float(global_js_map), int(global_js_map_rank), float(cons_js_map), int(cons_js_map_rank)])
    return (table_row, cons_js_score_array, residue_name_array, cons_js_rank_array)
