#!/usr/bin/env python

from itertools import product
import numpy as np
from os.path import expanduser
import pandas as pd
import sys
from typing import List, Tuple

from cx_analysis.skeleton import Skeleton


#####################
# dataframe_org.py
#
# Methods to extract and save summary data from a Connectome

def assemble_linkdf(C) -> pd.DataFrame:
    """
    link_df contains a row for each synaptic contact made between neurons in the Connectome
    :param C: Connectome
    :return link_df: DataFrame
    """
    df_rows = []
    skel_data = C.skel_data

    for pre_id, pre_sk in skel_data.items():
        assert (type(pre_sk) is Skeleton)
        out_links = pre_sk.out_links  # list containing a Dict for each synaptic link
        for l in out_links:
            post_id = l.get('post_skel')
            post_sk = C.skel_data.get(post_id, None)

            if post_sk is None:  # unidentified neurites (aka fragments)
                post_name = str(post_id)
                post_type = 'UNKNOWN'
                post_om = 'UNKNOWN'
            else:
                post_name = post_sk.name
                post_type = post_sk.subtype
                post_om = post_sk.group
            # TODO pd.Category this?
            df_rows.append({'pre_neuron': pre_sk.name,
                            'pre_type': pre_sk.subtype,
                            'pre_om': pre_sk.group,
                            'pre_skel': pre_id,
                            'post_neuron': post_name,
                            'post_type': post_type,
                            'post_om': post_om,
                            'post_skel': post_id,
                            'link_id': l.get('link_id'),
                            'cx_id': l.get('cx_id')})

    df = pd.DataFrame(data=df_rows, columns=['link_id', 'cx_id',
                                             'pre_neuron', 'pre_om', 'pre_type', 'pre_skel',
                                             'post_neuron', 'post_om', 'post_type', 'post_skel'])
    return df


def assemble_cxdf(C, linkdf) -> Tuple[pd.DataFrame, List, List]:
    """
    Longform DataFrame that has a row for each group of neurons/each connection type 
    
    requires link_df
    :param C: Connectome
    :return cxdf, inter, unknowns:
    """
    cx_types = [f"{pre}->{post}"
                for pre, post in product(C.cfg['subtypes'], C.cfg['subtypes'])]

    om_list = sorted([str(k) for k in C.grouping.keys()])

    counts = np.zeros((len(om_list), len(cx_types)), dtype=int)
    inter = []
    unknowns = []

    for ind, row in linkdf.iterrows():
        this_pre, this_post = (row['pre_type'], row['post_type'])
        if this_pre.upper() == 'UNKNOWN' or this_post.upper() == 'UNKNOWN':
            unknowns.append(row)
        elif row['pre_om'] != row['post_om']:
            inter.append(row)
        else:
            j = cx_types.index(f"{this_pre}->{this_post}")
            i = om_list.index(row['pre_om'])
            counts[i, j] += 1

    om_order = np.array([[om] * len(cx_types) for om in om_list]).reshape(-1)
    cx_order = np.tile(cx_types, len(om_list))
    print(f"om_order: {om_order.shape}, cx_order: {cx_order.shape}")
    pre_order, post_order = np.array([cx.split("->") for cx in cx_order]).T
    print(f"pre_order: {pre_order.shape}, post_order: {post_order.shape}")

    df = pd.DataFrame({'om': pd.Categorical(om_order),
                       'cx_type': pd.Categorical(cx_order),
                       'pre_type': pd.Categorical(pre_order),
                       'post_type': pd.Categorical(post_order),
                       'n_connect': np.ravel(counts)})
    df.loc[df['n_connect'] < 0, 'n_connect'] = np.nan

    return df, inter, unknowns


def extract_connector_table(linkdf: pd.DataFrame) -> pd.DataFrame:
    """
    Extract synaptic terminals from link dataframe
    TODO: split into two dataframes (one of summary, the other for cx's partner subtype breakdown)
    :param link_df:
    :return:
    """
    cx_data = dict.fromkeys(np.unique(linkdf['cx_id']))
    p_counts = dict.fromkeys(np.unique(linkdf['cx_id']))
    subtypes = sorted(np.unique(linkdf['post_type']))

    for cx, links in linkdf.groupby('cx_id'):
        cx_data[cx] = dict()
        # if something goes wrong here, li
        # presynaptic info
        cx_data[cx].update({'pre_om': np.unique(links['pre_om'])[0],
                         'pre_type': np.unique(links['pre_type'])[0],
                         'pre_neuron': np.unique(links['pre_neuron'])[0]})
        # Number of partners belonging to each subtype
        type_freq = links['post_type'].value_counts().to_dict()
        cx_data[cx].update({str(s): type_freq.get(s, 0) for s in subtypes})

        # Ommatidia of post partners (should be the same)
        partner_oms = np.unique(links['post_om'])
        partner_oms = partner_oms[partner_oms != 'UNKNOWN']
        if len(partner_oms) != 1:
            # raise Warning(f"Post-partners for connector {cx} belong to more than one ommatidia: {partner_oms}")
            cx_data[cx].update({'post_om': 'multi'})
        else:
            cx_data[cx].update({'post_om': str(partner_oms[0])})

    return pd.DataFrame(cx_data).T


def assemble_cxvectors(linkdf: pd.DataFrame, external: bool=True, excl_unknowns: bool=True) -> pd.DataFrame:
    """
    assemble_cxvectors
    Get a cxvectors dataframe from linkdf. 
    :param link_df: pd.DataFrame, longform containing a row for each synaptic contact, 
    :param external: bool, when True: synapses between neurons from different ommatidia are counted as a seperate category. 
    The postsynaptic (recipient) subtype for interommatidial connections will start with an 'e' 
    :param exclude_unknowns: bool, when True (default): synapses with an unidentified partner will not be counted. 
    When False: 'UNKNOWN' is treated like another subtype, giving the dataframe columns like 'R2R5->UNKNOWN' 
    
    
    Note: filtering certain connection types (e.g. discard connections with mean < 1.0) SHOULD NOT be done here
          (maybe even connection types with an UNKNOWN partner should be included in the output) 
    """
    # filter out connections to unidentified fragment
    if excl_unknowns:
        linkdf = linkdf.loc[((linkdf['pre_om'] != 'UNKNOWN') & (linkdf['post_om'] != 'UNKNOWN'))]
    
    oms = np.unique(linkdf['pre_om'])
    subtypes = np.unique([*linkdf['pre_type'], *linkdf['post_type']])
    ctypes = [f'{pre}->{post}' for pre, post in [p for p in product(subtypes, subtypes)]]
    # initialize df with all counts = 0
    df = pd.DataFrame(np.zeros((len(oms), len(ctypes)), dtype=int), index=oms, columns=ctypes)
    
    # if you want external connections to be listed in row of the pre neuron's home ommatidium, use
    # groupby 'pre_om' (R1R4->eL4). If you want them listed in the post neuron's home ommatidium, use 
    # groupby 'post_om' (eR1R4->L4)
    for om, rows in linkdf.groupby('pre_om'):
        for i, link in rows.iterrows():
            
            pre_type = link['pre_type']
            post_type = link['post_type']
            # local connections
            if link['pre_om'] == link['post_om']:
                df.loc[om, f'{pre_type}->{post_type}'] += 1
            # external connections 
            elif external:  # if pre_om != post_om
                if f'{pre_type}->e{post_type}' not in df.columns:
                    # make new ctype column when this external connection is first seen
                    df[f'{pre_type}->e{post_type}'] = np.zeros(len(oms), dtype=int)
                df.loc[om, f'{pre_type}->e{post_type}'] += 1
            else:
                continue # skip external connections if asked to
                
    return df    


def assemble_rhabdomere_df(full_df: pd.DataFrame) -> pd.DataFrame:
    # The excel file contains a table for each ommatidium
    data = []
    # TODO: extract other data
    om = []
    cell = []
    c5_offset = 0
    for i in range(0, 29):  # ommatidia
        # excel has 16 rows for each ommatidia, but only 13 have data
        this_range = full_df.loc[i * 16 + c5_offset: i * 16 + 13 + c5_offset].reset_index(drop=True)
        this_om = this_range.iloc[0, 0]
        if this_om == 'C5':  # C5 has an extra row 
            this_range = full_df.loc[i * 16: i * 16 + 14].reset_index(drop=True)
            rows = 11
            c5_offset += 1
        else:
            rows = 10
        # these cells contain the z-index where each measurement was taken
        z_st_cols = this_range.iloc[3, 13:]  

        assert(len(this_om) == 2)  # check om name
        assert(len(z_st_cols) == 9)  # check that there are 9 photoreceptors

        for ii, this_st in this_range.iloc[0, 4:13].items():  # subtypes 
            if '(' in this_st:  # some have the old subtype nomenclature in ()
                this_st = this_st.split('(')[0]
                
            this_st = this_st.strip().upper()
                
            z_col = z_st_cols[z_st_cols == this_st.lower()].index[0]
            
            # Add negative sign to stack index that start proximal 
            if this_range.iloc[2, 4] < this_range.iloc[3, 4]: # start/end index for R1
                z_inds = this_range.iloc[4:, z_col]
            else:
                z_inds = this_range.iloc[4:, z_col] * -1.0
            
#             if this_st == "R7'":
#                 this_st = 'R7p'
            #plus_offset = this_range.iloc[4:, ii] + offsets.loc[this_om, 1]
            data.append(pd.DataFrame({'om': [this_om]*rows, 
                                      'subtype': [this_st]*rows, 
                                      'z-index': z_inds, 
                                      'angle': this_range.iloc[4:, ii]}))

    df = pd.concat(data, ignore_index=True)
    df = df.astype({'om': str, 
                    'subtype': str,
                    'z-index': float, 
                    'angle': float})
    
    #compute some secondary results
    return rhabdomere_computations(df)

    
def rhabdomere_computations(df: pd.DataFrame) -> pd.DataFrame:
    
    from math import isnan
    
    for i, om_rows in df.groupby('om'):
        for ii, rows in om_rows.groupby('subtype'):
            n = 0
            for iii, row in rows.sort_values('z-index').iterrows():
                df.loc[iii, 'n'] = int(n)
                if n == 0:  # first measurement has diff of 0
                    df.loc[iii, '_diff'] = np.nan
                    df.loc[iii, 'diff'] = np.nan
                    df.loc[iii, 'cumulative'] = 0.0
                    df.loc[iii, '_angle'] = df.loc[iii, 'angle']
                    df.loc[iii, 'interval_len'] = np.nan # 0.0
                    df.loc[iii, 'interval_z'] = np.nan
                    previous = (iii, df.loc[iii, 'angle'])
                elif isnan(df.loc[iii, 'angle']):
                    #print(f"NaN found for {i}_{ii} measurement: {n}")
                    df.loc[iii, '_diff'] = np.nan
                    df.loc[iii, 'diff'] = np.nan
                    df.loc[iii, 'cumulative'] = np.nan
                    df.loc[iii, '_angle'] = np.nan
                    df.loc[iii, 'interval_len'] = np.nan
                    df.loc[iii, 'interval_z'] = np.nan
                    # previous = measurement before the NaN
                else:
                    df.loc[iii, '_diff'] = (df.loc[iii, 'angle'] - previous[1]) % 360.0
                    df.loc[iii, 'diff'] = df.loc[iii, '_diff'] - 360.0 * (df.loc[iii, '_diff'] > 180.0)
                    df.loc[iii, 'cumulative'] = df.loc[previous[0], 'cumulative'] + df.loc[iii, 'diff']
                    df.loc[iii, '_angle'] = df.loc[previous[0], '_angle'] + df.loc[iii, 'diff']
                    df.loc[iii, 'interval_len'] = (df.loc[iii, 'z-index'] - 
                                                   df.loc[previous[0], 'z-index']) * 8.0 / 1000.0
                    df.loc[iii, 'interval_z'] = df.loc[iii, 'z-index'] - df.loc[previous[0], 'z-index']
                    # because 1 px = 8/1000 microns
                    previous = (iii, df.loc[iii, 'angle'])
                n += 1
    df['diff_per_micron'] = df['diff'] / df['interval_len']
    df['cosine_sq'] = np.cos(df['angle']*(np.pi/180.0)) ** 2
    df['CCW_angle'] = [(360.0 + x) if x < 0 else x for x in df['angle']]

    return df
    
                        
            
    