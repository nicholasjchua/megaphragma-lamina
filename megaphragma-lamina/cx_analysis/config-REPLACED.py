#!/usr/bin/env python

from collections.abc import Mapping
import json
import os
import os.path
import logging
from glob import glob
from typing import List, Tuple, Dict
from typing import Dict
from collections import namedtuple

_Config = namedtuple(
    'Config',
    [
        'source',
        'cm_url',
        'cm_token',
        'p_id',
        'annot',
        'subtypes',
        'expected_n',
        'groupby',
        'annotator_initials',
        'min_cx',
        'save',
        'log',
        'out_dir',
        'restrict',
    ]
)

##################
# The Config object holds information representing experimental configuration,
# common to most code used in this library.

class Config(_Config):

    def cm_access(self) -> Tuple:
        return self.cm_url, self.cm_token, self.p_id


##################
# Class Methods

def parse_cfg_file(path: str="") -> Config:
    """
    parse_config_file
    :param path: str (optional), path to a json containing analysis parameters
    :return cfg: Dict
    """

    if path == "":  # the default file will cause it to crash rn
        fp = 'config/default_cfg.json'
    elif path.split('.')[-1] == 'json':
        fp = os.path.join(os.path.expanduser(path))
        path = os.path.split(fp)[0]
    else:  # if path to directory given
        fp = glob(os.path.join(os.path.expanduser(path), '*.json'))[0]
    print(fp)

    with open(fp) as cfg_file:
        cfg = json.load(cfg_file)
        
    print(path)

    # Check that the config file has all the necessary info
    ### URL TO CATMAID SERVER ###
    if cfg['cm_url'] == "":
        cm_url = os.environ['CM_URL']
    else:
        cm_url = cfg['cm_url']
        print("WARNING: Using Catmaid URL from config file. Remember to gitignore this if private.")
        
    ### USER ACCESS TOKEN FOR CATMAID ###
    if cfg['cm_token'] == "":
        cm_token = os.environ['CM_TOKEN']
    else:
        cm_token = cfg['cm_token']
        print("WARNING: User access token found in config file. Remember to gitignore this if private.")
    
    ### ANNOTATION ASSOCIATED WITH THE NEURONS YOU WANT TO QUERY ###
    if type(cfg['annot']) != str:
        raise Exception("A valid Catmaid annotation is required to fetch desired skeletons")
    else:
        annot = cfg['annot']
    
    ### CATMAID PROJECT ID ###
    if type(cfg['p_id']) != int:
        raise Exception("Project ID is an integer")
    else:
        p_id = cfg["p_id"]

    ### CX COUNT THRESHOLD (NOT USED, TODO: REMOVE) ###
    if type(cfg['min_cx']) != int:
        raise Exception("Connection count threshold is an integer")
    else:
        min_cx = cfg['min_cx']
    
    ### SAVE PREPROCESSED DATA ? ###
    if cfg['save'] is False:
        save = False
        out_dir = ""
        raise Warning("Data will not be saved based on options listed in config file")
    else:
        save = True
        # If no output dir given, save in location of cfg file
        if cfg['out_dir'] == "":
            out_dir = path
        else:
            out_dir = ""
    
    ### LOG (NOT YET IMPLEMENTED) ###
    if cfg['log'] is False:
        log = False
        raise Warning("Log will not be kept based on options listed in config file")
    else:
        log = True
    
    ### ANNOTATIONS USED TO CATEGORIZE NEURONS ###
    # Prespecified subtypes
    if type(cfg["subtypes"]) != list:
        raise Exception("Subtype categories need to be strings inside a list")
    else:
        subtypes = cfg["subtypes"]
    
    ### EXPECTED NUMBER OF EACH CATEGORY ###
    if len(subtypes) != len(cfg["expected_n"]):
        raise Exception("expected_n should be a the same length as subtypes")
    else:
        expected_n = cfg['expected_n']

    ### ANNOTATIONS USED TO CATEGORIZE NEURONS ###
    # ommatidium coord
    if cfg['groupby'] == 'om':
        groupby = 'om'
        annotator_initials = []  # don't need this if grouping by ommatidia
    elif cfg['groupby'] == 'annotator':
        if len(cfg['annotator_initials']) > 1:
            groupby = 'annotator'
            annotator_initials = cfg['annotator_initials']  # for validation experiment
        else:
            raise Exception("Annotator initials missing from config file (required when groupby=annotator)")
    else:
        raise Exception("Invalid 'groupby' argument. Needs to be 'annotator' or 'om'. Former also requires"
                        "a list of annotator initials in the config file")
    
    ### RESTRICT ANALYSIS ON A TAGGED SEGMENT OF THE SKELETON ###
    # Currently excludes all nodes/connections between the root node, and 'lamina_end'
    restrict = cfg['restrict_skeletons']
    # if cfg['restrict_skeletons']['restrict_tags'] == []:
    #     restrict =
    # else:
    #     restrict = cfg['restrict_skeletons']

    return Config(
        source=fp,
        cm_url=cm_url,
        cm_token=cm_token,
        p_id=p_id,
        annot=annot,
        subtypes=subtypes,
        expected_n=expected_n,
        groupby=groupby,
        annotator_initials=annotator_initials,
        min_cx=min_cx,
        save=save,
        log=log,
        out_dir=out_dir,
        restrict=restrict
    )

