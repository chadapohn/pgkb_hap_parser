#!/usr/bin/env python

import glob
import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [atoi(c) for c in re.split(r"(\d+)", text)]

def get_sort_allele_definition_table_file(allele_definition_path):
    sort_allele_definition_table_file = []
    for allele_definition in glob.glob(allele_definition_path + "/*.xlsx"):
        sort_allele_definition_table_file.append(allele_definition)
    sort_allele_definition_table_file.sort(key=natural_keys)
    return sort_allele_definition_table_file
