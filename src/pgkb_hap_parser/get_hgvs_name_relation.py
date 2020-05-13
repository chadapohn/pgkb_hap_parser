#!/usr/bin/env python

def get_hgvs_name_relation(hgvs, haplotype_cell):
    row = haplotype_cell.shape[0]
    col = haplotype_cell.shape[1]
    hgvs_name_relation = {}
    for i in range(row):
        if 0 != i:
            name_relation = []
            for j in range(col):
                if 0 != j:
                    if "nan" != str(haplotype_cell.iloc[i, j]):
                        name_relation.append(haplotype_cell.iloc[0, j])
            hgvs_name_relation[hgvs[i - 1]] = name_relation
    return hgvs_name_relation
