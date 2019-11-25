#!/usr/bin/env python
from os import path
import pandas as pd
import re

HAPLOTYPE_TABLE_DIR = path.join(path.dirname(__file__), "haplotype-tables")
GENE_ROW = 0
GENE_COL = 0
GENE_PATTERN = r"^GENE:\s*(\w+)$"

CHROM_ROW = 3
CHROM_COL = 1
CHROM_PATTERN = r"^.*N\w_\s*(\d+)\.\d+"

def parse_gene(gene_cell):
    gene = None
    match = re.match(GENE_PATTERN, gene_cell)
    if match:
        gene = match.group(1)
    return gene

def parse_chrom(chrom_cell):
    chrom_num = None
    match = re.match(CHROM_PATTERN, chrom_cell)
    if match:
        chrom_num = int(match.group(1))
        if chrom_num == 23:
            chrom_name = "chrX"
        elif chrom_num == 24:
            chrom_name = "chrY"
        else:
            chrom_name = "chr" + str(chrom_num)
    return chrom_num
        
if __name__ == "__main__":
    definition_file = path.join(HAPLOTYPE_TABLE_DIR, "G6PD_allele_definition_table.xlsx")
    definition_table = pd.read_excel(definition_file, header=None)
    gene = parse_gene(definition_table.iloc[GENE_ROW][GENE_COL])
    chrom = parse_chrom(definition_table.iloc[CHROM_ROW][CHROM_COL])