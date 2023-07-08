#!/usr/bin/env python

from typing import Tuple

def hex_to_rgb(hex: str) -> Tuple[int, int, int]:
    # FIXME: This looks to be another copy of the same function defined elsewhere in this library, returning a
    # 3-member rather than 4-member tuple
    h = hex.lstrip('#')
    assert(len(h) == 6)
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

