#!/usr/bin/env python

import re

def get_rsid(rsid_cell):
    rsid = []
    for cell in rsid_cell.values.tolist():
        match_rsid = re.match(r"^rs\d+$", str(cell))
        if match_rsid:
            rsid.append(match_rsid.group())
        else:
            rsid.append("")
    return rsid
