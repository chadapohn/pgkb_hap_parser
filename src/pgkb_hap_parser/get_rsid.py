#!/usr/bin/env python

import re

def get_rsid(rsid_cell):
    rsid = []
    for cell in rsid_cell.values.tolist():
        match_rsid = re.findall(r"rs\d+", str(cell))
        if match_rsid:
            if len(match_rsid) > 1:
                rsid.append(match_rsid)
            else:
                rsid.append(match_rsid[0])
        else:
            rsid.append("")
    return rsid
