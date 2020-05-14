#!/usr/bin/env python

def get_name_relation_to_hgvs(hgvs, haplotype_cell):
    row = haplotype_cell.shape[0]
    col = haplotype_cell.shape[1]
    name_relation_to_hgvs = {}
    for i in range(row):
        if 0 != i:
            relation_to_hgvs = []
            for j in range(col):
                if 0 != j:
                    if "nan" != str(haplotype_cell.iloc[i, j]):
                        relation_to_hgvs.append(hgvs[j - 1])
            name_relation_to_hgvs[haplotype_cell.iloc[i, 0]] = relation_to_hgvs
    return name_relation_to_hgvs
