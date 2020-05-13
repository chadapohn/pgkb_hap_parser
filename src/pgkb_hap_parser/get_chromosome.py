#!/usr/bin/env python

import re

def get_chromosome(chromosome_cell):
    match_chromosome = re.search(r"chromosome\s(\w+)", str(chromosome_cell))
    if match_chromosome:
        return "chr" + match_chromosome.group(1)
    else:
        return
