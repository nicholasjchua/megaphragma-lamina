#!/usr/bin/env python

import numpy as np
from typing import Dict, Union, Tuple, List
from scipy.spatial import distance

from cx_analysis.catmaid_queries import get_root_id, node_with_tag, skel_compact_detail, node_coords
#from cx_analysis.config import Config
import sys
from os.path import expanduser


def nodes_betwixt(skel_id: str, cfg, restrict_tags: Union[str, Tuple], node_data: List=None,
                       invert: bool=False) -> Union[List[str], Tuple]:
    """
    Get a list of node_id that are between two tags in the skeleton. Invert will return the complement, the 
    nodes outside the tagged segment
    TODO: allow this to take node_ids for start and end instead
    :param skel_id: Skeleton ID
    :param restrict_tags: str or tuple of two strings. Giving just one will define the segment as root -> <tag>
    :param nodes: (optional) list of node IDs so they aren't fetched again
    :param invert: If true, returns the nodes OUTSIDE the tagged segment.
    :return:
    """
    root_id = get_root_id(skel_id, cfg)

    if node_data is None:
        node_data = skel_compact_detail(skel_id, cfg)

    if type(restrict_tags) is str:
        start = root_id
        end = node_with_tag(skel_id, root_id, restrict_tags, cfg)  # str
        
    elif len(restrict_tags) == 1:
        start = root_id
        end = node_with_tag(skel_id, root_id, restrict_tags[0], cfg)
    elif len(restrict_tags) == 2:
        start = node_with_tag(skel_id, root_id, restrict_tags[0], cfg)
        end = node_with_tag(skel_id, root_id, restrict_tags[1], cfg)
    else:
        raise Exception("More than two restrict_tags given")

    dist = check_dist(start, end, cfg)
    nodes_within = traverse_nodes(node_data, int(start), int(end))
    
    # returns only the node_ids as a list of strs 
    if invert: 
        return [str(n[0]) for n in node_data if n[0] not in nodes_within]
    else:
        return [str(n) for n in nodes_within]


def traverse_nodes(node_data: List, start: int, end: int) -> List:
    """
    Recursive walk through node_list from 'start' node, collecting node IDs in a list.
    Function will split when a branch is hit (i.e. node has two children)
    Ends when when an arbor terminates (has no children)
    :param node_list: List, of nodes to traverse
    :param start:
    :param end:
    :return node_ids: A List of node IDs  
    """
    current = start
    node_ids = []
    
    if current != end:
        node_ids.append(current)
        children = [n[0] for n in node_data if n[1] == int(current)]
        if children == []:  # end of branch
            return node_ids
        else:
            for c in children:
                deeper_nodes = traverse_nodes(node_data, start=c, end=end)
                node_ids.extend(deeper_nodes)

    return node_ids # as list of ints, but converted to strs in the wrapper above 


def segment_skeleton(skel_id: str, cfg, node_data: List=None, restrict_tags: Union[Tuple, str]=None, restrict_nodes: List=None, 
                     verbose: bool=False) -> Dict:
    """
    Get a dict describing the branch structure of a skeleton:
    Each key is the ID of a branch (first node after a fork). This key points to a list of nodes between this branch node and a fork or end
    List associated with each branch starts with the parent of the branch node
    Calls _segment_skeleton() to do the recursive branch traversal
    :param skel_id: ID of skeleton to be segmented
    :param cfg: Config file for Catmaid access 
    :node_data: Passing a node data 2D array will avoid called skel_compact_detail
    :restrict_tags: Tags used to exclude certain nodes from the segmentation with node_betwixt(invert=True) if restrict_nodes are NOT passed
    :restrict_nodes: Nodes excluded from segmentation (e.g. nodes outside the lamina)
    :verbose: If true, reports on the way nodes are restricted, and gives a summary of results
    """
    root_id = int(get_root_id(skel_id, cfg))

    if node_data is None:
        node_data = skel_compact_detail(skel_id, cfg)
    
    if (restrict_nodes is not None) and (restrict_tags is None):
        if verbose: 
            print(f"Skeleton {skel_id}: {len(restrict_nodes)} nodes excluded")
        restrict_nodes = [int(r) for r in restrict_nodes]
        node_data = [n for n in node_data if int(n[0]) not in restrict_nodes]
    elif restrict_tags is not None:
        restrict_nodes = nodes_betwixt(skel_id, cfg, restrict_tags, node_data, invert=True)
        restrict_nodes = [int(r) for r in restrict_nodes]
        node_data = [n for n in node_data if int(n[0]) not in restrict_nodes]
    else:
        node_data = node_data
    
    seg_map = {root_id: []}
    return _segment_skeleton(node_data, segments=seg_map, current=root_id, parent_branch=root_id)


def _segment_skeleton(node_data: List, segments: Dict, current: int, parent_branch: int) -> Dict:
    """
    Recursive step of segment_skeleton()
    - Base case: the current node has no children (is a leaf). Current is added to its parent branch's list in the dict. Dict is returned. 
    - If current has > 1 child, it is a branch point. Add current to the parent branch's list. 
      A new key and list is added to 'branches'. 
    - If current has 1 child, append ID to the list in {parent_branch: []}, continue to its children without changing the current 'parent_branch' 
    """
    
    #these_branches = branches
    #b = branches
    children = [n[0] for n in node_data if n[1] == int(current)]
    #b.setdefault(last_branch, []).append(current)
    
    if len(children) == 0:  # End of segment
        segments[parent_branch].append(int(current))
        return segments
    
    elif len(children) > 1:  # Branch point
        segments[int(parent_branch)].append(current) 
        #segments.update({int(current): [int(current)]})
        for child in children:
            # new dict entry for each child, current added to start of their lists
            segments.update({child: [current]})
            deeper_segs = _segment_skeleton(node_data, current=int(child), 
                                         segments=segments, parent_branch=int(child))
            segments.update(deeper_segs) 
        return segments
    
    elif len(children) == 1: # 1 child, continue on, parent remains the same
        segments[parent_branch].append(current)
        deeper_segs = _segment_skeleton(node_data, current=int(children[0]), 
                                        segments=segments, parent_branch=int(parent_branch))
        segments.update(deeper_segs) 
        return segments
    
    else:
        raise Exception(f"Current node {current} has no children in node_list")
    
    
def dist_two_nodes(n1: str, n2: str, cfg=None, coord_map: Dict=None) -> float:
    
    if (coord_map is None) & (cfg is not None):
        # Will do a catmaid query
        coord1 = np.array(node_coords(n1, cfg), dtype=float)
        coord2 = np.array(node_coords(n2, cfg), dtype=float)
    elif coord_map is not None:
        coord1 = np.array(coord_map[n1], dtype=float)
        coord2 = np.array(coord_map[n2], dtype=float)
    else:
        raise Exception("cfg is needed if node_coords are not given")
    
    dist = distance.euclidean(coord1, coord2)
    return dist


def check_dist(n1: str, n2: str, cfg) -> float:
    """
    Raises an exception if two nodes are too close (i.e. lamina_end and root)
    :param n1: str, node_id of root
    :param n2:
    :param cfg:
    :return dist: float, if sufficiently far, returns the distance in nm
    """
    dist = dist_two_nodes(n1, n2, cfg)
    if np.abs(dist) < 1000.0:
        raise Exception(f'Nodes too close ({dist: .01f}nm)')
    else:
        return dist



def find_central_segment(skel_id: str, end_tag: str, cfg, node_data: List=None, restrict_nodes: List=None, verbose: bool=False) -> List:
    """
    Get a list of nodes between a tagged ending and the root 
    """
    root_node = int(get_root_id(skel_id, cfg))
    target_node = int(node_with_tag(skel_id, root_node, end_tag, cfg))
    
    print(target_node)
    
    if node_data is None:
        node_data = fetch_node_data(skel_id, cfg)
    if restrict_nodes is not None:
        node_data = [n for n in node_data if n[0] not in restrict_nodes] 
        
    parent_map = dict()
    for row in node_data:
        parent_map.update({row[0]: row[1]})  # node: node's parent
    
    central_segment = []
    current = target_node
    # Collect intervening nodes from target to root
    while current != root_node:
        central_segment.append(current)
        current = parent_map[current]  # parent of current becomes the new current
        
    if verbose:
        print(f"skel_id {skel_id} - central segment: {len(central_segment)} nodes, total: {len(node_list)} nodes")
    
    return central_segment



def find_end_points(node_list: List) -> List:
    """
    Get IDs of end nodes (leaf nodes) of the skeleton
    """
    parent_nodes = np.array(node_list).T[1]
    end_ids = [n[0] for n in node_list if n[0] not in parent_nodes]
    
    return end_ids


def measure_path_lengths(branch_list: Dict, cfg=None, node_data: List=None) -> Dict:
    """
    Measure the path length of each branch segment in a skeleton
    """
    if node_data is not None:
        coord_map = {row[0]: row[2:5] for row in node_data}
    elif (node_data is None) & (cfg is None):
        raise Exception("cfg is needed if node_list is not passed")
    else:
        coord_map = None
    # loop across branch segments   
    results = dict.fromkeys(branch_list.keys())
    for b, seg in branch_list.items():
        results[b] = seg_length(seg, cfg, coord_map)
        
    return results

def seg_length(seg: List, cfg=None, coord_map: Dict=None) -> float:
    """
    Path length of a list of connected nodes
    """
    this_len = 0.0
    for n1, n2 in zip(seg[0 : -1], seg[1: ]):
        this_len += dist_two_nodes(n1, n2, cfg, coord_map)
    
    return this_len 


def measure_seg_distances(branch_list: Dict, cfg=None, node_data: List=None) -> Dict:
    """
    Measure the euclidean distance between the first and last node of each branch segment in a 
    skeleton
    """
    if node_data is not None:
        coord_map = {data[0]: data[3:6] for data in node_data}
    elif (node_data is None) & (cfg is None):
        raise Exception("cfg is needed if node_list is not passed")
    else:
        coord_map = None
        
    results = dict.fromkeys(branch_list.keys())
    for b, seg in branch_list.items():
        results[b] = seg_distance(seg, cfg, coord_map)
    
    return results

def seg_distance(seg: List, cfg=None, coord_map: Dict=None) -> float:
    """
    Euclidean distance of the first and last points of a segment
    """
    return dist_two_nodes(seg[0], seg[-1], cfg, coord_map)

