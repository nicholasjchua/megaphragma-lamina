from typing import Dict, Tuple, Union, List, Iterable
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import networkx as nx
from pprint import pprint

def hexplot(node_data: Union[Dict, pd.DataFrame], decimals: int=0, 
            var_lim: Iterable=None, c: Iterable=None, lc: str='k',
            ax: plt.Axes=None, scale_factor: float=0.1):#float=0.015):
    """
    Plot a map of the hexagonal lattice of the megaphragma compound eye. Each ommatidium will be labelled and coloured
    according to the strings and (r, g, b, a) values passed with node_data. This will also take a 1-col dataframe indexed by om
    TODO: fix dataframe input options, type of cmap as an arg, use plt.scatter instead of all the networkx functions? 
    :param node_data: Dict, {om: {'label': str,
                                 {'outline': matplotlib line spec, e.g. '-',
                                 {'colour':  (rgba)}}
    :param edge_data:
    :param ax: plt.Axes, if None, will plot to current Axis
    :param scale_factor: float, controls the spacing between nodes
    :param n_rgba: Default node colour (if node_data: colours is None)
    :param e_colour: Default edge colour (if edge_data is used)
    :return:
    """
    if c is None: # default color
        c = (0.2, 0.2, 0.2, 1)
        
    if not isinstance(node_data, Dict):
        node_data = __from_series(node_data, c=c, var_lim=var_lim, decimals=decimals)
    
    G, pos = generate_lattice()
    nx.set_node_attributes(G, pos, name='pos')
    # Handle labels/data
    node_colours = []
    node_outline = []
    node_labels = {}
    name_to_ind = {}
    for nx_ind, data in G.nodes(data=True):
        this_om = hex_to_om(data['pos'])
        name_to_ind.update({this_om: tuple(nx_ind)})
        nd = node_data.get(this_om, {})
        label = nd.get('label', this_om)
        
        node_colours.append(nd.get('colour', c))
        node_outline.append(nd.get('outline', '-'))
        
        node_labels.update({nx_ind: nd.get('label', this_om)})

    if ax is None:
        ax = plt.gca()

    pos = scale_distance_between(pos, scale_factor)

    nx.draw(G, pos, alpha=1.0, node_color=node_colours, node_size=300, #scale_factor * 18000,
            node_shape='H', ax=ax)
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=5, font_color=lc, ax=ax)
#     ax.set_xmargin(0)
#     ax.set_ymargin(0)
    
    xs, ys = ([v[0] for k, v in pos.items()], 
              [v[1] for k, v in pos.items()])
    
    ax.set_xlim([min(xs) - scale_factor, max(xs) + scale_factor])
    ax.set_ylim([min(ys) - scale_factor, max(ys) + scale_factor])
    ax.set_aspect('equal')
    
    return ax


def hexplot_TEST(node_data, decimals: int=0, 
                 var_lim: Iterable=None, c: object=None, lc: str='k',
                 #edge_data: Dict=None, edge_c='r',  # EDGE DATA NOT IMPLEMENTED YET
                 ax: plt.Axes=None, scale_factor: float=0.015):
    """
    Plot a map of the hexagonal lattice of the megaphragma compound eye. Each ommatidium will be labelled and coloured
    according to the strings and (r, g, b, a) values passed with node_data. This will also take a 1-col dataframe indexed by om
    TODO: fix dataframe input options, type of cmap as an arg, use plt.scatter instead of all the networkx functions? 
    :param node_data: Dict, {om: {'label': str,
                                 {'outline': matplotlib line spec, e.g. '-',
                                 {'colour':  (rgba)}}
    :param edge_data:
    :param ax: plt.Axes, if None, will plot to current Axis
    :param scale_factor: float, controls the spacing between nodes
    :c: Default node colour (if node_data: colours is None)
    :param e_colour: Default edge colour (if edge_data is used)
    :return:
    """
    if c == None: # default color
        c = (0.2, 0.2, 0.2, 1)
        
    if not isinstance(node_data, Dict):
        node_data = __from_series(node_data, c=c, var_lim=var_lim, decimals=decimals)
        
    if ax == None:
        ax = plt.gca()
        
    om_list = sorted([str(om) for om in node_data.keys()])
    pos = [om_to_hex(o) for o in om_list] # 2D figure coords of each om
    node_colours = []#dict.fromkeys(om_list)
    node_outline = []#dict.fromkeys(om_list)
    node_labels = []#dict.fromkeys(om_list)
        
    #name_to_ind = 
    for om, xy in zip(om_list, pos):
        
        if node_data[om].get('label', None) == None:
            label = om
        elif isinstance(node_data.get('label'), (int, float)):
            label = str(round(node_data.get('label'), decimals))
        else:
            label = node_data.get('label')
                        
        if (node_data[om].get('colour') == None):
            fill_c = c
        else:
            fill_c = node_data[om].get('colour')
            
        x, y = (xy[0] * 0.01, xy[1] * 0.01)
        #y = xy[1] * 0.01
        
        ax.scatter(xy[0], xy[1], marker='H', color=fill_c, s=100)
        ax.annotate(label, xy, fontsize=8, color='w', ha='center', va='center')
        
        
    #ax.set_xlim((-30, 4))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    #plt.axis('off')
    ax.set_aspect('equal')
    
    return ax




def generate_lattice() -> Tuple:
    """
    hex_lattice
    :param spacing: float, factor that scales the distance between points in pos
    :return: networkx Graph, with nodes representing each ommatidia in *Megaphragma*
    :return: pos, Dict of 2D coordinates for plotting (if scaling != 1.0, these will be different from the coordinates
             in G)
    """
    G = nx.generators.triangular_lattice_graph(4, 10)
    G.remove_node((0, 4))  # delete the most ventro-medial
    pos = nx.get_node_attributes(G, name='pos')  # a dict of node_ref: fig_coordinates

    # Rotate 270 so that vertical axis is anatomical ventral -> dorsal, horizontal is medial -> lateral
    for om in pos.keys():
        tmp = pos[om]
        pos[om] = ((tmp[1] * -1), tmp[0])
        pos[om] = (pos[om][0] + float(4.0 * (np.sqrt(3.0) / 2.0)), pos[om][1] + 0.5)
        if om in [(0, 2), (0, 3), (1, 4)]:  # some nodes need to be repositioned further
            pos[om] = (pos[om][0] + float(np.sqrt(3) / 2), pos[om][1] - 0.5)
    edges = list(G.edges())
    for e in edges:
        G.remove_edge(e[0], e[1])
    #G = nx.set_node_attributes(G, pos, name='pos')

    return G, pos


def __from_series(X: pd.Series, c: Iterable, var_lim: Iterable=None, center_cmap_at: str=None, decimals: int=0) -> Dict:
    """
    __from_series
    This function is called when hexplot() receives a pandas Series instead of 
    a dict for node_data. Right now, the series has to be indexed by 'om'. 
    The different options for the color transfer functions need to be tested. This is 
    currently good if you don't need a cmap centered on a certain value, i.e. min -> max, white -> red
    1. 
    """
    
    from matplotlib.colors import LinearSegmentedColormap
    
    if isinstance(X, pd.Series):
        var_name = X.name
    else:
        raise Exception('node_data is not a Pandas Series')
        
        
    ### MAKE COLORMAP ACCORDING TO OPTIONS AND TRANSFORM DATA ###
    if center_cmap_at is None: # cmap goes from white to c (min val -> max val)
        cm = LinearSegmentedColormap.from_list(name='mycm', colors=[(1.0, 1.0, 1.0), c])
        if var_lim is None:  # range uses max and min of series
            trans_vals = (X - X.min()) / (X.max() - X.min())
        else:  # custom range
            trans_vals = (X - var_lim[0]) / (var_lim[1] - var_lim[0])
    ### NOT TESTED ###
    elif center_cmap_at == 'mean':
        if len(c) == 3 and isinstance(c[0], (int, float)):
            cm = LinearSegmentedColormap.from_list(name='mycm', colors=[c, (1.0, 1.0, 1.0), c])
            raise Warning("Pass 2 RGB tuples for a cmap that diverges into two colors for min and max. " +
                          "Currently only one RGB value given")
        elif isinstance(c[0], Iterable) and len(c[0]) == 3: 
            cm = LinearSegmentedColormap.from_list(name='mycm', colors=[c[0], (1.0, 1.0, 1.0), c[1]])
        else:
            raise Exception("c argument invalid, should be a tuple containing RGB value, " + 
                            "or a list of 2 RBG tuples for a divergent colormap.")
        trans_vals = abs(X/X.mean()) - 0.5
    
#     elif 'mid' in center_cmap_at:
#         mid_point = (X.min() + X.max()) * 0.5
#         if len(c) == 3 and isinstance(c[0], (int, float)):
#             cm = LinearSegmentedColormap.from_list(name='mycm', colors=[c, (1.0, 1.0, 1.0), c])
#             raise Warning("Pass 2 RGB tuples for a cmap that diverges into two colors for min and max. " +
#                           "Currently only one RGB value given")
#         elif isinstance(c[0], Iterable) and len(c[0]) == 3: 
#             cm = LinearSegmentedColormap.from_list(name='mycm', colors=[c[0], (1.0, 1.0, 1.0), c[1]])
#         else:
#             raise Exception("c argument invalid, should be a tuple containing RGB value, " + 
#                             "or a list of 2 RBG tuples for a divergent colormap.")
#         trans_vals = (X/mid_point) + 0.5
    ##################
    else:
        raise Exception('Something went wrong')
    
    node_data = dict()
    for om, val in X.items():
        this_node = {'colour': cm(trans_vals.loc[om])}
        if isinstance(val, int):  # integer
            this_node.update({'label': f"{val: .0f}"})
        else: # non-int
            v = np.round(val, decimals)
            this_node.update({'label': f"{val}"})
        
        node_data[om] = this_node
            
    return node_data


def hex_to_om(xy: Iterable) -> str:
    """
    Convert a set of figure coordinates to a letter-digit ommatidia ID
    :param position: 2-tuple of floats indicating the ommatidia's figure coordinates
    :return col_row: str, e.g. 'B2'
    """
    assert(len(xy) == 2)
    col_num = np.rint(4 - (2.0/np.sqrt(3.0) * xy[0]))  # A is most anterior, evaluates to 0
    row_num = int(np.floor(xy[1]) + col_num/2.0)
    col_letter = chr(int(col_num) + 65)
    return str(col_letter) + str(row_num)


def om_to_hex(om: str, scale_factor: float=0.1) -> Tuple:
    """
    Convert letter-digit ommatidia ID (e.g. 'A4') to figure coordinates
    :param om: length-2 string describing the ommatidia's coordinate on the eye
    :return x, y: x and y coordinates for 2D figures
    """
    assert(len(om) == 2)
    col_num = float(ord(om[0]) - 65.0)
    x = (-1.0 * np.sqrt(3.0) / 2.0 * col_num + 4.0) * scale_factor
    y = (float(om[1]) - (0.5 * col_num)) * scale_factor

    return x, y


def scale_distance_between(positions, scalar):
    """
    Takes a dict of figure coordinates and scales the distance between them
    :param positions: dict(node_ref: (x, y))
    :param scalar: distances between nodes is multiplied by this number
    :return: positions: the results of the scaling operation
    """
    for this_node, this_pos in positions.items():
        positions[this_node] = (this_pos[0] * scalar, this_pos[1] * scalar)
    return positions

