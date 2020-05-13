#!/usr/bin/env python

import re

def get_gene(gene_cell):
    match_gene = re.match(r"^GENE:\s*(\w+)\s*$", str(gene_cell))
    if match_gene:
        return match_gene.group(1)
    else:
        return
