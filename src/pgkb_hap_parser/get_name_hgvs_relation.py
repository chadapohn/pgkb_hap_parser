#!/usr/bin/env python

def get_name_hgvs_relation(hgvs, haplotype_cell):
    row = haplotype_cell.shape[0]
    col = haplotype_cell.shape[1]
    name_hgvs_relation = {}
    for i in range(row):
        if 0 != i:
            hgvs_relation = []
            for j in range(col):
                if 0 != j:
                    if "nan" != str(haplotype_cell.iloc[i, j]):
                        hgvs_relation.append(hgvs[j - 1])
            name_hgvs_relation[haplotype_cell.iloc[i, 0]] = hgvs_relation
    return name_hgvs_relation
