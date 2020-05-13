#!/usr/bin/env python

import re
import copy

iupac = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]",
    "-": "del"
}

def get_haplotype_name_variant_extract_iupac(haplotype_cell):
    haplotype = haplotype_cell.values.tolist()
    # get haplotype nan to reference
    for i in range(len(haplotype)):
        for j in range(len(haplotype[i])):
            if "nan" == str(haplotype[i][j]):
                haplotype[i][j] = haplotype[0][j]
    # get haplotype transform iupac
    for i in range(len(haplotype)):
        for j in range(len(haplotype[i])):
            if haplotype[i][j] in iupac:
                haplotype[i][j] = iupac.get(haplotype[i][j])
    # get haplotype extract iupac
    haplotype_extract_iupac = []
    for i in range(len(haplotype)):
        haplotype_extract_iupac.append(haplotype[i])
        for j in range(len(haplotype[i])):
            match_iupac = re.match(r"^\[(\D+)\]$", str(haplotype[i][j]))
            if match_iupac:
                del haplotype_extract_iupac[-1]
                variant_iupac = match_iupac.group(1)
                variant_iupac = list(variant_iupac)
                for k in range(len(variant_iupac)):
                    haplotype_iupac = copy.deepcopy(haplotype)
                    haplotype_iupac[i][j] = variant_iupac[k]
                    haplotype_extract_iupac.append(haplotype_iupac[i])
    # get name_extract_iupac, variant_extract_iupac
    name_extract_iupac = []
    variant_extract_iupac = []
    for i in range(len(haplotype_extract_iupac)):
        name_extract_iupac.append(haplotype_extract_iupac[i][0])
        variant_extract_iupac.append(haplotype_extract_iupac[i][1:])
    return haplotype_extract_iupac, name_extract_iupac, variant_extract_iupac
