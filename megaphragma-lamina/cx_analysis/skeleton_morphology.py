#!/usr/bin/env python

from typing import List, Dict, Tuple, Union
import numpy as np
import pandas as pd
import json

from cx_analysis.node_ops import segment_skeleton, find_central_segment, measure_path_lengths, measure_seg_distances
from cx_analysis.connectome import Connectome
from cx_analysis.catmaid_queries import get_root_id

"""
skeleton_morphology.py
Methods to compute and save morphology data. 
:Date: 13-Apr-2021
:Authors: Nicholas Chua 
"""

def run_morphology_analysis(C: Connectome, skel_ids: List[str], 
                            restrict_tags: Union[Tuple, str]=None, 
                            save_file: str=None, verbose: bool=False) -> Tuple[Dict, Dict, Dict, Dict, Dict]:
    """
    Run a series of morphology measurements for specified skeletons in a Connectome object

    Parameters
    ----------
    :param C: Connectome, object defined in cx_analysis.connectome
    :param skel_ids: List[str], a list of skeleton IDs (str) to run morphology analysis
    :param restrict_tags: Union[Tuple, str], CATMAID tags pointing to start and end nodes to 
    constrain analysis. When one str is given, uses the root as the start point 
    (has only been tested for cases with one restrict node, 'lamina_end')
    :param save_path: str, optional, full file path to save json of data.
    :param verbose: bool, print information on each skeleton for debugging

    :return: tuple, containing four dictionaries: 
        segments,
        central_segs,
        seg_lengths,
        seg_distances
    """
    segments = dict()
    central_segs = dict()
    seg_lengths = dict()
    seg_distances = dict()
    strahler = dict()

    for s in skel_ids:  # TODO make it loop through C instead 
        s = str(s) # skel_id keys in C.skel_data are str
        if s not in list(C.skel_data.keys()):
            raise Exception(f"Skeleton with ID {s} not found in Connectome.skel_data")
        else:
            data = C.skel_data[s]
            # this returns a dict where each entry is a list of contiguous node IDs
            segments[s] = segment_skeleton(s, cfg=C.cfg, 
                                           node_data=data.skel_nodes, 
                                           restrict_nodes=data.r_nodes, 
                                           verbose=verbose)
            # List of nodes along the 'backbone'. This will break if end_tag is not str
            central_segs[s] = find_central_segment(s, end_tag=restrict_tags,
                                                   cfg=C.cfg, 
                                                   restrict_nodes=data.r_nodes, 
                                                   node_data=data.skel_nodes)
            # measure path lengths of each segment
            seg_lengths[s] = measure_path_lengths(segments[s], 
                                                  cfg=C.cfg, 
                                                  node_data=data.skel_nodes)
            # measure distance between the start and end of each segment
            seg_distances[s] = measure_seg_distances(segments[s], 
                                                     cfg=C.cfg, 
                                                     node_data=data.skel_nodes)

            strahler[s] = strahler_order(segments[s], data.skel_nodes,
                                         data.r_nodes)
        if verbose:
            print(f"Morphology data computed for {s}")

    # Save results as json
    if save_file is not None:
        #print(save_file)
        save_morphology_data(save_file, segments, central_segs, seg_lengths, seg_distances, strahler)

    return segments, central_segs, seg_lengths, seg_distances, strahler


def save_morphology_data(fn: str, segments: Dict, central_segs: Dict,   # kwags when finalized
                         seg_lengths: Dict, seg_distances: Dict, strahler: Dict) -> None:
    """
    Save results of run_morphology_analysis as a json
    """
    if fn[-5:] != '.json':
        fn = fn + '.json'

    results = {'segments': segments, 
                'central_segs': central_segs, 
                'seg_lengths': seg_lengths, 
                'seg_distances': seg_distances, 
                'strahler': strahler}

    with open(fn, 'x') as fh:
        json.dump(results, fh)
    print(f"Morphology data saved as {fn}")


def strahler_order(segments: Dict, node_data: List, r_nodes: List=[]) -> Dict:
    """
    Compute the strahler order of each branch point in a skel segments dict
    """
    if r_nodes:
        # the ids saved in C.skel_data.r_nodes is str
        r_nodes = [int(n) for n in r_nodes]
        node_data = [n for n in node_data if n[0] not in r_nodes]

    #root = find_root_node(node_data, skel_data, cfg)
    root = find_root_node(node_data)
    leaves = find_leaf_nodes(node_data) 

    # precompute map of each branch point's next upstream or downstream bps
    pbps = find_parent_branches(segments)
    cbps = find_child_branches(segments)

    # keep track of the points where we need to calculate SO 
    # leaf nodes cannot also be parent branch points
    #points = list(set(list(leaves + list(pbps.keys()) + [root])))
    points = list(set(list(leaves + list(cbps.keys()))))
    n_points = len(points)

    so = dict()
    so_calculated = []
    # Deal with leaf nodes first, SO = 1
    for leaf in leaves:
        if leaf == root:
            continue
        else:
            so[leaf] = 1
            so_calculated.append(leaf)
            # Append this leaf's SO to a list of each its parent
            so.setdefault(pbps[leaf], []).append(1)

    while len(so_calculated) < (n_points - 1): 
        for p in points:

            if not so.get(p, None): # bp has no entry in SO,
                continue  # hold off until SO computed for children
            elif type(so[p]) is int:  # SO already calculated, move on
                continue
            # This point has a list with 
            elif len(cbps[p]) == len(so[p]):
                so[p] = max(so[p]) + 1
                so_calculated.append(p)
                # Update this point's parent (unless root, which has no parent)
                if p == root: 
                    continue
                else:
                    so.setdefault(pbps[p], []).append(so[p])
            else:
                continue

    if type(so[root]) is list:
        so[root] = max(so[root]) + 1

    return so


# def _strahl(segments: Dict, current_bp: int, node_data: List,  
#             leaves: List, results: Dict):

#     child_nodes = [n[0] for n in node_data if n[1] == int(current_bp)]
#     child_bps = [segments[c] for c in children]
#     child_orders = []
#     for cbp in child_bps:
#         if cbp in leaves:
#             child_orders.append(1)
#         else:
#             #child_results = _strahl(

#     results[current_segment] = max([_strahl(c) for c in child_branches])


# def reverse_segments(segments: Dict) -> Dict:
#     """
#     Reverse the dictionaries of segmented nodes so that each segment is indexed
#     by the downstream branch point or terminal node (centripetal format)
#     """
#     #print({int(seg[-1]): [int(s) for s in seg.reverse()] for k, seg in segments.values()})
#     return {int(seg[-1]): [int(s) for s in seg.reverse()] for k, seg in segments.items()}


def find_parent_branches(segments: Dict) -> Dict:
    """
    Map each branch point or terminal node to its next parent branch
    INPUT THE ORIGINAL SEG DICT, not the reversed one
    """
    #r_seg = reverse_segments(segments)
    #return {r[k]: v[-1] for k, v in r_seg.items()}
    pbs = {int(seg[-1]): int(seg[0]) for k, seg in segments.items()}
    #print(pbs)
    return pbs


def find_child_branches(segments: Dict) -> Dict:
    #print(segments)
    child_bps = dict()
    for seg in segments.values():
        child_bps.setdefault(seg[0], []).append(seg[-1])
    return child_bps

def find_root_node(node_data: List) -> int:
    """
    Quick and easy way to find the root node from node_data
    doesn't require making a server query like catmaid_queries.get_root_id()
    :param node_data: 2D list, each line is a node, 
                      e.g. [2156, None, 4, 36192.0, 73024.0, 25376.0, 0.0, 5]
                           [node_id, parent_id, ?, x, y, z, ?, ?]
    :return: int, node_id of root
    """

    root = [int(n[0]) for n in node_data if n[1] is None] # no parent
    if len(root) != 1:  
        raise Exception(f"Found {len(root)} root nodes in node_list")
    else:
        return int(root[0])

def find_leaf_nodes(node_data: List) -> List[int]:
    """
    Quick and easy way to find the leaf nodes in node_data
    doesn't require making a server query
    :param node_data: 2D list, each line is a node, 
                      e.g. [2156, None, 4, 36192.0, 73024.0, 25376.0, 0.0, 5]
                           [node_id, parent_id, ?, x, y, z, ?, ?]
    :return: List[int], node_id of root
    """
    parent_nodes = np.array(node_data).T[1]
    return [int(n[0]) for n in node_data if n[0] not in parent_nodes]


def precomputed_children(node_data: List) -> Dict:
                
    the_map = dict()
    for this_node in np.array(node_data).T[0]:
        the_map[this_node] = [n[0] for n in node_data if n[1] == int(this_node)]

    return the_map
