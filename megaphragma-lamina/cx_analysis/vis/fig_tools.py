from typing import Union, Dict, Tuple, List
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def hex_to_rgb(hex: str) -> Tuple:
    h = hex.lstrip('#')
    assert(len(h) == 6)

    tup = tuple(int(h[i:i+2], 16) / 256.0 for i in (0, 2, 4))
    return tup


def linear_cmap(n_vals: int, max_colour: Union[Tuple, str], min_colour: Union[Tuple, str] = (1.0, 1.0, 1.0),
                mp_colour: Union[Tuple, str]=None) -> LinearSegmentedColormap:
    if mp_colour is None:
        clist = np.array([min_colour, max_colour])
    else:
        clist = np.array([min_colour, mp_colour, max_colour])

    cm = LinearSegmentedColormap.from_list('mycm', clist, N=n_vals)
    return cm


def anatomical_axes(x: float, y: float, arrow_len: float, x_text: str="anterior", y_text: str="dorsal") -> None:
    """
    TODO: Use ax.annotate instead so that a predefined Axes object can be passed into the func
    :param x:
    :param y:
    :param arrow_len:
    :param x_text:
    :param y_text:
    :return None:
    """
    plt.arrow(x, y, arrow_len, 0.0, fc='k', ec='k', head_width=0.03, head_length=0.05)
    plt.text(x + arrow_len + 0.1, y, x_text, fontsize=12)
    plt.arrow(x, y, 0.0, arrow_len, fc='k', ec='k', head_width=0.03, head_length=0.05)
    plt.text(x, y + arrow_len + 0.1, y_text, fontsize=12)

    
def midpoint_2d(point1: Union[List, Tuple], point2: Union[List, Tuple]) -> Tuple:
    """
    Get the midpoint of two points in 2D space
    :param point1:
    :param point2:
    :return mp: tuple of 2D coordinates
    """
    return (point1[0] + point2[0])/2.0, (point1[1] + point2[1])/2.0

def shortened_st(labels: List) -> List:
    """
    Convert each label in a list of subtype labels to their shortened form. 
    R2R5 etc are unchanged. 'LMC_X' becomes LX. 'centri' becomes 'Am'. R7p becomes "R7'" 
    Warning: when shortened names are used, the single quote mark in R7' might cause errors
    """
    
    
    short_labels = []
    for l in labels:
        if l == 'R7p':
            short_labels.append("R7'")
        elif l[0:3] == 'LMC':
            short_labels.append(l[0] + l[-1])
        elif l == 'centri':
            short_labels.append('Am')
        elif l == 'R_duet':
            short_labels.append('R2+R5')
        elif l == 'R_quartet':
            short_labels.append('R1+R3+R4+R6')
        else:
            short_labels.append(l)
            
    return short_labels
            

def subtype_cm(fp: str="../cx_analysis/vis/lamina_palette.json", rgb: bool=False, include_spr=True) -> Dict:
    """
    Get map of 'subtype' to 'hex colour' from json file
    :param fp: str, File containing colormap
    :param rgb: bool, Returns a Tuple of (R,G,B) instead of hexadecimal
    :param include_spr: Include individual short photoreceptors as keys, e.g. 'R1' in addition to 'R1R4'
    :return palette: Dict, with either hex or tuples as values
    """
    with open(fp) as f:
        p = json.load(f)
    
    if include_spr:
        p.update({'R1': p['R1R4'], 'R4': p['R1R4'],
                        'R2': p['R2R5'], 'R5': p['R2R5'],
                        'R3': p['R3R6'], 'R6': p['R3R6']})
    if rgb:
        return {k: hex_to_rgb(v) for k, v in p.items()}
    else:
        return p
    
    