#!/usr/bin/env python
from os import path
import pandas as pd
import re

HAPLOTYPE_TABLE_DIR = path.join(path.dirname(__file__), "haplotype-tables")
GENE_ROW = 0
GENE_COL = 0
GENE_PATTERN = r'^GENE:\s*(\w+)$'

def parse_gene(gene_cell):
    match = re.match(GENE_PATTERN, gene_cell)
    gene = None
    if match:
        gene = match.group(1)
    return gene
        
if __name__ == "__main__":
    definition_file = path.join(HAPLOTYPE_TABLE_DIR, "G6PD_allele_definition_table.xlsx")
    definition_table = pd.read_excel(definition_file, header=None)
    gene = parse_gene(definition_table.iloc[GENE_ROW][GENE_COL])
    