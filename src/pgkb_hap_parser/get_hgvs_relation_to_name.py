#!/usr/bin/env python

def get_hgvs_relation_to_name(hgvs, haplotype_cell):
    haplotype_cell = haplotype_cell.T
    row = haplotype_cell.shape[0]
    col = haplotype_cell.shape[1]
    hgvs_relation_to_name = {}
    for i in range(row):
        if 0 != i:
            relation_to_name = []
            for j in range(col):
                if 0 != j:
                    if "nan" != str(haplotype_cell.iloc[i, j]):
                        relation_to_name.append(haplotype_cell.iloc[0, j])
            hgvs_relation_to_name[hgvs[i - 1]] = relation_to_name
    return hgvs_relation_to_name
