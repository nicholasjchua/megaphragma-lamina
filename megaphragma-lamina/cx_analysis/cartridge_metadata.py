#!/usr/bin/env python

from typing import Dict, List, Tuple, Union
import matplotlib.pyplot as plt
from datetime import date


def lamina_subtypes(rgb: bool=False) -> Dict:
    """
    Colour for each neuron category, parsed from the correlation matrix labels
    :param labels: List, with labels of connection types in the format: 'pre->post'
    :param p: str, use the 'pre' or 'post' category
    :return:
    """
    c_code = {'R1R4': '#A2D16E',
              'R3R6': '#dcff82',
              'R2R5': '#587932',
              'R7': '#843f8a',
              'R7p': '#6000E9',
              'R8': '#b603fc',
              'LMC_1': '#EA33F7',
              'LMC_2': '#8D2E4F',
              'LMC_3': '#EB512E',
              'LMC_4': '#F3AF3D',
              'LMC_N': '#e08ba6',
              'centri': '#0307fc'}
    if rgb:
        return {k: hex_to_rgb(v) for k, v in c_code.items()}
    else:
        return c_code


def hex_to_rgb(hex: str) -> Tuple[float, float, float, float]:
    """ Converts a hex string like #AAACCC" to a tuple of the 3 colour bits composing it, plus a trailing value of 1 (default alpha) """
    h = hex.lstrip('#') + 'ff'
    return tuple(float(int(h[i: i+2], 16))/255.0 for i in (0, 2, 4, 6))



def annotator_id() -> Dict[str, str]:
    """ Lookup table for who worked on what bodies in our project """
    # TODO: Remove this before making the repo public
    initials = {'A2': 'MT',
               'A3': 'NC',
               'A4': 'BC',
               'A5': 'SV',
               'B0': 'BC',
               'B2': 'SV',
               'B3': 'SV',
               'B4': 'MT',
               'B5': 'BC',
               'B6': 'SV',
               'C3': 'SV',
               'C4': 'NC',
               'C5': 'MT',
               'C6': 'SV',
               'D2': 'SV',
               'D3': 'SV',
               'D4': 'BC',
               'D5': 'BC',
               'D6': 'MT',
               'D7': 'SV',
               'E4': 'BC',
               'E5': 'BC',
               'E6': 'SV',
               'E7': 'SV'}

    return initials


def r7_inputs_handcount() -> Dict[str, int]:
    """
    Data on 24 ommatidia from Sep 5 2019
    :return:
    """
    # TODO: Remove this before making the repo public
    r7_inputs = {'A2': 1,
               'A3':  9,
               'A4': 18,
               'A5': 24,
               'B0': 3,
               'B2': 0,
               'B3': 4,
               'B4': 4,
               'B5': 20,
               'B6': 19,
               'C3': 2,
               'C4': 7,
               'C5': 12,
               'C6': 24,
               'D2': 2,
               'D3': 0,
               'D4': 3,
               'D5': 8,
               'D6': 9,
               'D7': 18,
               'E4': 0,
               'E5': 7,
               'E6': 16,
               'E7': 17}

    return r7_inputs

def ret_clusters() -> Dict[str, List[str]]:
    """
    DRA cluster based on connectivity data, VRA are R8vol > R7vol
    :return ret_clust: Dict mapping each om to (dorsal, central, or ventral)
    """
    # TODO: Remove this before making the repo public
    ret_clust = {'dra': ['A4', 'A5', 'B5', 'B6', 'C5', 'C6',
                         'D6', 'D7', 'E6', 'E7'],
                 'vra': ['A0', 'A1', 'A2', 'B0', 'B3', 'C3',
                         'D2', 'D3', 'E4'],
                 'v_trio': ['B1', 'C1', 'C2']}
    return ret_clust

def anastasia_clusters() -> Dict[str, str]:
    """
    Based on Anastasia's clustering of ommatidia according to rhabdom measurements (190723)
    """
    # TODO: Remove this before making the repo public
    r7_inputs = {'A2': 'classic',
               'A3':  'classic',
               'A4': 'dra minor',
               'A5': 'dra',
               'B0': 'ventral',
               #'B1'and 'A0': 'ventral',
               'B2': 'classic',
               'B3': 'classic',
               'B4': 'classic',
               'B5': 'dra',
               'B6': 'dra',
               'C3': 'classic',
               'C4': 'classic',
               'C5': 'dra minor',
               'C6': 'dra',
               'D2': 'classic',
               'D3': 'classic',
               'D4': 'classic',
               'D5': 'classic',
               'D6': 'dra',
               'D7': 'dra',
               'E4': 'classic',
               'E5': 'classic',
               'E6': 'dra minor',
               'E7': 'dra'}

    return r7_inputs


def anastasia_clust_cmap() -> Dict:
    # TODO: Remove this before making the repo public

    clusters = anastasia_clusters()

    key = {'ventral': 'r',
           'classic': 'orange',
           'dra minor': 'violet',
           'dra': 'darkviolet'}

    node_data = {k: {'label': k,
                     'colour': key[v]} for k, v in clusters.items()}

    return node_data

def annotator_cmap() -> Dict:
    # TODO: Remove this before making the repo public

    initials = annotator_id()

    rgb = {'MT': 'r',
           'NC': 'c',
           'SV': 'y',
           'BC': 'g'}

    node_data = {k: {'label': k,
                     'colour': rgb[v]} for k, v in initials.items()}

    return node_data


def om_days_since() -> Dict:
    # TODO: Remove this before making the repo public

    om_dates = {'A2': '12-2018',
                'A3': '10-2018',
                'A4': '07-2019',
                'A5': '06-2019',
                'B0': '06-2019',
                'B2': '01-2019',
                'B3': '10-2018',
                'B4': '01-2019',
                'B5': '08-2019',
                'B6': '03-2019',
                'C3': '12-2018',
                'C4': '10-2018',
                'C5': '10-2018',
                'C6': '10-2018',
                'D2': '08-2019',
                'D3': '03-2019',
                'D4': '01-2019',
                'D5': '03-2019',
                'D6': '06-2019',
                'D7': '08-2019',
                'E4': '06-2019',
                'E5': '07-2019',
                'E6': '07-2019',
                'E7': '08-2019'}

    now = date.today()
    days = {k: (now - date(int(v[3:7]), int(v[0:2]), 1)).days
            for k, v in om_dates.items()}

    return days


def om_days_cmap() -> Dict:
    # TODO: Remove this before making the repo public

    days = om_days_since()

    node_data = {k: {'id': k,
                     'label': v}
                 for k, v in days.items()}

    max_months = max([v['label'] for k, v in node_data.items()])
    cmap = plt.get_cmap('Greens', max_months)

    for k, v in node_data.items():
        node_data[k]['colour'] = cmap(v['label'])
        if v['label'] > 252:
            node_data[k]['label'] = str(v['label']) + '*'

    return node_data


def r7_inputs_cmap() -> Dict:
    # TODO: Remove this before making the repo public

    counts = r7_inputs_handcount()
    cmap = plt.get_cmap('RdPu', max([round(v) for k, v in counts.items()]))
    node_data = {k: {'colour': cmap(int(v)),
                     'id': k,
                     'label': str(round(v))} for k, v in counts.items()}


    return node_data
