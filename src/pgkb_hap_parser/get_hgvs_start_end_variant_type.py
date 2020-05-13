#!/usr/bin/env python

import re

def get_hgvs_start_end_variant_type(hgvs_cell):
    hgvs = hgvs_cell.values.tolist()
    start = []
    end = []
    variant_type = []
    for cell in hgvs:
        match_snp = re.match(r"^g\.(\d+)([A-Z]>([A-Z]|[A-Z](\/[A-Z])*))*$", str(cell))
        match_ins = re.match(r"^g\.(\d+)_?(\d*)ins.+$", str(cell))
        match_del = re.match(r"^g\.(\d+)_?(\d*)del.+$", str(cell))
        if match_snp:
            start.append(match_snp.group(1))
            end.append(match_snp.group(1))
            variant_type.append("SNP")
        elif match_ins:
            start.append(match_ins.group(1))
            end.append(match_ins.group(2))
            variant_type.append("INS")
        elif match_del:
            start.append(match_del.group(1))
            if match_del.group(1) and not match_del.group(2):
                end.append(match_del.group(1))
            elif match_del.group(1) and match_del.group(2):
                end.append(match_del.group(2))
            else:
                end.append("")
            variant_type.append("DEL")
        else:
            start.append("")
            end.append("")
            variant_type.append("")
    return hgvs, start, end, variant_type
