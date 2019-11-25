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

VARIANT_ROW = 3
VARIANT_COL = 2
IN_OR_DEL_PATTERN = r"^[cgp]\.(\d+\_\d+).*$"
SNP_PATTERN = r"^[cgp]\.(\d+).*$"

RSID_ROW = 5
RSID_COL = 2
RSID_PATTERN = r"rs\d+"

def parse_gene(gene_cell):
    gene = None
    match = re.match(GENE_PATTERN, gene_cell)
    if match:
        gene = match.group(1)
    return gene

def parse_chrom(chrom_cell):
    chrom_name = None
    match = re.match(CHROM_PATTERN, chrom_cell)
    if match:
        chrom_num = int(match.group(1))
        if chrom_num == 23:
            chrom_name = "chrX"
        elif chrom_num == 24:
            chrom_name = "chrY"
        else:
            chrom_name = "chr" + str(chrom_num)
    return chrom_name

def parse_variants(variant_cells, num_variants):
    starts = []
    ends = []
    chrom_hgvs_names = []
    for var_idx in range(VARIANT_COL, num_variants):
        chrom_hgvs_name = variant_cells.loc[var_idx]
        match_in_or_del_pattern = re.match(IN_OR_DEL_PATTERN, chrom_hgvs_name)
        match_snp_pattern = re.match(SNP_PATTERN, chrom_hgvs_name)
        if match_in_or_del_pattern:
            start = int(match_in_or_del_pattern.group(1).split("_")[0]) - 1
            end = int(match_in_or_del_pattern.group(1).split("_")[1])
        elif match_snp_pattern:
            end = int(match_snp_pattern.group(1))
            start = end - 1
        else:
            chrom_hgvs_name = None
            start = None
            end = None
        chrom_hgvs_names.append(chrom_hgvs_name)
        starts.append(start)
        ends.append(end)

    return chrom_hgvs_names, starts, ends

def parse_rsids(rsid_cells, num_variants):
    rsids = []
    for rsid_idx in range(RSID_COL, num_variants):
        rsid = rsid_cells.loc[rsid_idx]
        if pd.notna(rsid):
            match = re.match(RSID_PATTERN, rsid)
            if match:
                rsid = match.group()
        else:
            rsid = None
        rsids.append(rsid)

    return rsids

if __name__ == "__main__":
    definition_file = path.join(HAPLOTYPE_TABLE_DIR, "G6PD_allele_definition_table.xlsx")
    definition_table = pd.read_excel(definition_file, header=None)
    gene = parse_gene(definition_table.iloc[GENE_ROW][GENE_COL])
    chrom = parse_chrom(definition_table.iloc[CHROM_ROW][CHROM_COL])
    num_variants = definition_table.iloc[CHROM_ROW].count() - 1
    chrom_hgvs_names, starts, ends = parse_variants(definition_table.iloc[VARIANT_ROW, VARIANT_COL:num_variants+VARIANT_COL], num_variants)
    rsids = parse_rsids(definition_table.iloc[RSID_ROW, RSID_COL:num_variants+RSID_COL], num_variants)
    