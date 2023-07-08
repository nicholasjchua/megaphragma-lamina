#!/usr/bin/env python

from cx_analysis.catmaid_queries import *
#from cx_analysis.config import Config
from cx_analysis.node_ops import nodes_betwixt


class Skeleton:

    def __init__(self, skel_id: str, name: str, group: str, cfg):

        self.skel_id: str = skel_id
        self.name: str = name
        self.group: str = group
        self.cfg = cfg
        self.cartridge: str = name[2:4]  # which ommatidia this is from, from neuron_name, not 'ommatidia_XY' annotation
        self.subtype: str = self.__assign_subtype()  # check cfg to decide if categories are excl
        self.skel_nodes: List = skel_compact_detail(skel_id, cfg)

        if self.subtype in cfg['restrict_skeletons']['restrict_for']:
            print(f"{self.name}: RESTRICT")
            self.r_nodes: List = nodes_betwixt(skel_id, cfg, cfg['restrict_skeletons']['restrict_tags'], invert=True)
            print(f"{len(self.r_nodes)} nodes are in the restricted region")
        else:
            self.r_nodes: List = []
        self.out_cx, self.out_links, self.r_cx = cx_in_skel(skel_id, cfg, r_nodes=self.r_nodes)


    def __assign_subtype(self, allow_nulltype: bool=False) -> str:

        excl_categories = set(self.cfg['subtypes'])
        annotations_found = set(annot_in_skel(self.skel_id, self.cfg))

        intersect = list(excl_categories & annotations_found)
        if len(intersect) > 1:
            raise Exception(f"Skeleton {self.skel_id} subtype ambiguous: {intersect}")
        elif len(intersect) == 0:
            if allow_nulltype:
                return ''
            else:
                raise Exception(f"Skeleton {self.name} does not belong to any category in cat_list")
        else:
            return str(intersect[0])



