#!/usr/bin/env python

import pdb
from typing import List, Tuple, Union, Dict, Sequence
import requests
import sys
from os.path import expanduser

#from cx_analysis.config import Config

####################
# REQUEST WRAPPERS #
####################

def do_get(apipath: str, cfg) -> Tuple[bool, str]:
    """
    Wraps get requests. Performs specific action based on the operation specificed by apipath
    :param apipath: API path for a specific get operation, e.g. '/annotations/'
    :return response: True if the request was successful
    :return results: A json of the results if sucessful, a string if not.
    """
    print(type(cfg))
    print(cfg['cm_url'])
    p_url = cfg['cm_url']
    token = cfg['cm_token']
    p_id = cfg['p_id']
    
    path = p_url + "/" + str(p_id) + apipath
    result = requests.get(path, headers={'X-Authorization': 'Token ' + token})

    try:
        jresult = result.json()
    except ValueError:
        jresult = None

    if jresult is not None:
        if 'type' in jresult:  # API doesn't provide another way to get at a python-objects-structured parse
            if jresult['type'] == "Exception":
                print("exception info:")
                print(jresult['detail'])
                return False, "Something went wrong"
            
    if result.status_code == 200:
        if jresult is not None:
            return True, jresult
        else:
            raise Exception(f"Did not return json, but also not an error, text is {result.text}")
    else:
        return False, f"Something went wrong with {apipath}, return code was {result.status_code}"


def do_post(apipath: str, postdata: Dict, cfg) -> Tuple[bool, str]:
    """
    Wraps post requests. Performs specific action based on the operation specificed by apipath
    and the fields in postdata

    :return response: True if the request was successful
    :return results: A json of the results if successful, a string if not.
    """

    p_url = cfg['cm_url']
    token = cfg['cm_token']
    p_id = cfg['p_id']
    
    path = p_url + "/" + str(p_id) + apipath
    result = requests.post(path, data=postdata, headers={'X-Authorization': 'Token ' + token})

    try:
        jresult = result.json()
    except ValueError:
        jresult = None

    if jresult is not None:
        if 'type' in jresult:
            if jresult['type'] == "Exception":
                print("exception info:")
                print(jresult['detail'])
                return False, "Something went wrong"

    if result.status_code == 200:
        if jresult is not None:
            return True, jresult
        else:
            raise Exception(f"Did not return json, but also not an error, text is {result.text}")
    else:
        return False, f"Something went wrong with {apipath}, return code was {result.status_code}"


# ANNOTATION QUERIES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def skels_in_annot(annot: Union[int, str], cfg) -> Tuple[List[str], List[str]]:
    """
    Given an annotation ID, produce lists of skeleton IDs and neuron names
    :param: annot: annotation to query neurons
    :param:cfg: Config, object stores analysis configs
    :return: neuron_names: Lists of (skel_id, neuron_name)
    """

    if type(annot) is str:
        annot_id = annot_to_id(annot, cfg)
    else:
        annot_id = annot

    op_path = "/annotations/query-targets"
    post_data = {"annotated_with": annot_id,
                 "types": ["skeleton", "neuron"]}
    res_code, data = do_post(op_path, post_data, cfg)

    if len(data["entities"]) == 0:
        raise Exception(f"Entities annotated with annotation ID: {annot_id} not found")
    skeleton_ids = [str(entity.get("skeleton_ids")[0]) for entity in data["entities"]]
    neuron_names = [str(entity.get("name")) for entity in data["entities"]]

    if None in (skeleton_ids or neuron_names):
        raise Exception("Entities are missing fields")

    return skeleton_ids, neuron_names

def annot_in_skel(skel_id: str, cfg) -> List[str]:
    """
    Fetch list of annotations associated with the skeleton ID, raises exception if no annotations found
    :param cfg: Dict, of analysis options
    :param skel_ids: str,the numerical skeleton ID
    :return: annot_list, List, of annotations
    """
    op_path = "/annotations/forskeletons"
    post_data = {"skeleton_ids": skel_id}
    res_code, data = do_post(op_path, post_data, cfg)

    annot_list = list(data["annotations"].values())
    if annot_list is []:
        raise Exception(f"{skel_id} has no annotations")
    else:
        return annot_list

def annot_to_id(annot: str, cfg) -> int:
    """
    Search project for an annotation, get its numerical ID
    :param token: Catmaid API token
    :param p_id: Project ID
    :param annot: str, Annotation
    :returns annot_id: int
    """

    op_path = "/annotations/"
    res_code, data = do_get(op_path, cfg)
    data = data['annotations']
    annot_id = None
    for this_annot in data:
        if this_annot['name'] == annot:
            annot_id = this_annot['id']
    if annot_id is None:
        raise Exception(f"The annotation: {annot} does not exist in project: {cfg.p_id}")
    else:
        return annot_id


# TREENODE QUERIES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_root_id(skel_id: str, cfg, verbose: bool=False) -> str:
    """
    Get the node ID corresponding to the root of a skeleton_id
    :param skel_id: str
    :return node_id: str
    """
    op_path = f"/skeletons/{skel_id}/root"
    res_code, root = do_get(op_path, cfg)
    node_id = root.get('root_id', None)

    if node_id is None:
        raise Exception(f"Root node not found for skeleton: {skel_id}")
    else:
        if verbose:
            print(f"{skel_id} - root: {node_id}")
        return str(node_id)


def fetch_node_data(node_id: str, cfg) -> List:
    """
    Get data associated with a node ID
    get treenode/{node}/compact-detail gives node data in a slightly diff order than skeletons/{skel}/compact-detail
    :param node_id: str
    :param cfg:
    :return: List, of data corresponding to the node
    """
    op_path = f"/treenodes/{node_id}/compact-detail"
    res_code, node_data = do_get(op_path, cfg)
    # TODO what does this look like again?
    if type(node_data) is list:
        return node_data
    else:
        raise Exception(f"Could not find node with ID: {node_id}")

def node_with_tag(skel_id: str, root_id: str, tag_regex: str, cfg, first: bool=True, verbose: bool=False) -> Union[str, List]:
    """
    Returns the node_id of the first node in the skeleton tagged with 'tag_regex'

    Note: the api call returns a list of nodes in ascending distance from root_id, this function returns the one
    nearest to root. The tag could also be a regular expression.
    :param skel_id: Skeleton ID
    :param root_id: ID of root node. Used to sort tagged nodes by distance
    :param cfg: Config object
    :param tag_regex: Tag you want to query
    :return: node_id: str node ID of the tagged treenode. If first=False, returns a list of [nodeID,nodeData] in
    ascending order of distance from root.
    """
    # TODO Check the full output of this API call
    op_path = f"/skeletons/{skel_id}/find-labels"
    post_data = {"treenode_id": int(root_id),
                 "label_regex": str(tag_regex)}
    res_code, data = do_post(op_path, post_data, cfg)

    if len(data) == 0:
        raise Exception(f"Skeleton {skel_id} does not have a node tagged with {tag_regex}")
    elif first:
        if verbose:
            print(f"First node found in {skel_id} with {tag_regex} is {data[0][0]}")
        return str(data[0][0])
    else:
        if verbose:
            print(f"nodes with tag: {data}")
        return data


def node_coords(node_id: str, cfg) -> Tuple:
    """
    Get the x, y, z coordinates of a node using fetch_node_data,
    get treenode/{node}/compact-detail gives node data in a slightly diff order than skeletons/{skel}/compact-detail
    :param node_id:
    :return x, y, z
    """
    x, y, z = fetch_node_data(node_id, cfg)[2: 5]
#     if (x < 20) or (y < 20) or (z < 20):
#         raise Exception("Fetched data might be other node data, not coords. Check the result of fetch_node_data")
    return x, y, z


# SKELETON QUERIES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def skel_compact_detail(skel_id: str, cfg) -> List:
    """
    Returns a list of treenodes for a given skeleton
    The 'post' version of this API call doesn't seem to work, so I don't know how to have 'with_connectors' etc.
    nodes is a list of lists: [[nodes], [connectors], {nodeID: [tags]}], but currently the latter two are empty

    :param skel_id: Skeleton ID
    :param cfg: Config object
    :return: nodes[0], the first element which contains the node data
    """

    op_path = f"/skeletons/{skel_id}/compact-detail"

    res_code, nodes = do_get(op_path, cfg)
    return nodes[0]


# CONNECTOR QUERIES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def cx_in_box(x: int, y: int, z: int, l: int, cfg, v_size: int=8) -> Tuple[Dict, Dict]:
    """
    Get a Dict of connectors (pre syn) and postsynaptic nodes bounded by a
    cube starting at coordinates x, y, z with side length l (voxels).
    voxel size (default = 8) needs to be specified as the API takes 'real world' coordinates (nanometers)
    Returns voxel coordinates
    """

    op_path = f"/connectors/in-bounding-box"
    post_data = {'minx': x * v_size,
                 'miny': y * v_size,
                 'minz': z * v_size,
                 'maxx': (x + l) * v_size,
                 'maxy': (y + l) * v_size,
                 'maxz': (z + l) * v_size,
                 'with_locations': True,
                 'with_links': True}

    res_code, data = do_post(op_path, post_data, cfg)
    # API returns a list of lists. Each inner list contains the data for either a presynaptic or postsynaptic point.
    # field[0] = connector ID, field[1:4] = (x, y, z) in world coords, field[4] = skel_id of connector, field[-1] = 15 if post, 16 if pre 
    # Currently only interested in the coordinates and connectivity of pre/post nodes in the cube, skel_id can also be obtained
    # from this call if needed 

    if len(data) < 1 or data is None:
        raise Exception("Specified cube does not contain any synapses, make sure you are passing voxel coordinates") 
    else:
        pre_coords = dict()
        post_coords = dict()

        for cx in data:
            if cx[-1] == 15:  # presynaptic entry
                pre_coords.update({cx[0]: [float(c) / float(v_size) for c in cx[1:4]]})
            elif cx[-1] == 16:  # postsynaptic entry
                node_xyz = node_coords(cx[7], cfg)
                node_xyz = [float(c) / float(v_size) for c in node_xyz]
                post_coords.setdefault(cx[0], []).append(node_xyz)

        return pre_coords, post_coords

def cx_in_skel(skel_id: str, cfg, r_nodes: List) -> Tuple:
    """
    Get a Dict of all connectors PRESYNAPTICally associated with a neuron, and the associated link data
    """

    op_path = "/connectors/"
    post_data = {"skeleton_ids": skel_id,
                 "relation_type": "presynaptic_to", # doesn't seem to do anything (links w relation=15 still present)
                 "with_tags": True,
                 "with_partners": True}

    res_code, data = do_post(op_path, post_data, cfg)
    # data['partners'] is how you get the cx_id: [links] dictionary

    r_connectors = set()
    tmp_connector_data = dict()
    tmp_link_data = []
    for cx_id, p_data in data["partners"].items():
        for link in p_data:
            if link[3] == 15 and str(link[1]) in r_nodes:
                # print(f"Found a restricted connector {cx_id}")
                r_connectors.add(cx_id)
            elif link[3] == 16:
                tmp_link_data.append({'link_id': str(link[0]),
                                      'pre_skel': skel_id,
                                      'post_skel': str(link[2]),
                                      'post_node': str(link[1]),
                                      'cx_id': cx_id})
                tmp_connector_data.setdefault(cx_id, []).append(link)
    # Filter out restricted ones
    if len(r_connectors) > 0:
        link_data = [l for l in tmp_link_data if l['cx_id'] not in r_connectors]
        connector_data = {c: data for c, data in tmp_connector_data.items() if c not in r_connectors}
    else:
        link_data = tmp_link_data
        connector_data = tmp_connector_data

    return connector_data, link_data, list(r_connectors)

 
def cx_coords(cx_id: str, cfg) -> Tuple:
    """
    cx_coords
    Method to get the (x, y, z) coordinates associated with a connector_id
    """
    op_path = f"/connectors/{cx_id}"
    res_code, data = do_get(op_path, cfg)
    return data['x'], data['y'], data['z']

