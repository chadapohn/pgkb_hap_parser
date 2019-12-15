#!/usr/bin/env python
from os import path, listdir
import pandas as pd
import re
import numpy as np
from pprint import pprint
import json

HAPLOTYPE_TABLE_DIR = path.join(path.dirname(__file__), '..', "definition", "xlsx")
OUT_DIR = path.join(path.dirname(__file__), '..', "definition", "tsv") 

######
JSON_DIR = path.join(path.dirname(__file__), '..', "definition", "json") 

GENE_ROW = 0
GENE_COL = 0
GENE_PATTERN = "^GENE:\s*(\w+)$"

CHROM_ROW = 3
CHROM_COL = 0
CHROM_PATTERN = r"^.*N\w_\s*(\d+)\.\d+"

VARIANT_ROW = 3
VARIANT_COL = 1
SNP_PATTERN = r"^[cgp]\.(\d+)[A-Z]>[A-Z]"
SNP_INS_PATTERN = r"^[cgp]\.(\d+)_?(\d*)ins[A-Z]*"
SNP_DEL_PATTERN = r"^[cgp]\.(\d+)_?(\d*)del[A-Z]*"
SNP_REPEAT_PATTERN = r"^[cgp]\.(\d+)"
# SV_DEL_PATTERN = r"^delGene$"

RSID_ROW = 5
RSID_COL = 1
RSID_PATTERN = r"rs\d+"

HAP_ROW = 7
HAP_COL = 0


#######################
positions = []
alleles = None
hap_names = None

def parse_gene(gene_cell):
    gene_cell = gene_cell.strip()

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
    global positions
    starts = []
    ends = []
    chrom_hgvs_names = []
    var_types = []
    for var_idx in range(VARIANT_COL, num_variants + 1):
        chrom_hgvs_name = variant_cells.loc[var_idx]

        matched_snp = re.match(SNP_PATTERN, chrom_hgvs_name)
        matched_ins = re.match(SNP_INS_PATTERN, chrom_hgvs_name)
        matched_del = re.match(SNP_DEL_PATTERN, chrom_hgvs_name)
        matched_repeat = re.match(SNP_REPEAT_PATTERN, chrom_hgvs_name)

        # print(chrom_hgvs_name, matched_snp, matched_ins, matched_del)

        if matched_snp:
            position = matched_snp.group(1)

            ######
            positions.append(int(position))

            start = int(position) - 1
            end = int(position)
            var_types.append("SNP")
        elif matched_ins:
            position = matched_ins.group(1)
            if "_" in position:
                start, end = position.split("_")

                #####
                positions.append(int(start))

                start = int(start) - 1
                end = int(end)
            else:
                #####
                positions.append(int(position))

                start = int(position) - 1
                end = int(position)
            var_types.append("INS")
        elif matched_del:
            print(matched_del)
            position = matched_del.group(1)
            if "_" in position:
                start, end = position.split("_")

                #####
                positions.append(int(start))

                start = int(start) - 1
                end = int(end)
            else:
                #####
                positions.append(int(position))

                start = int(position) - 1
                end = int(position)
            var_types.append("DEL")
        elif matched_repeat:
            position = matched_repeat.group(1)

            #####
            positions.append(int(position))
            
            start = int(position) - 1
            end = int(position)
            var_types.append("SNP") # TODO: Temporary Fix
        else:
            #####
            positions.append("")
            chrom_hgvs_name = None
            start = None
            end = None
            var_types.append("")
        chrom_hgvs_names.append(chrom_hgvs_name)
        starts.append(start)
        ends.append(end)

    print(str(len(chrom_hgvs_names)), str(len(positions)))
    return chrom_hgvs_names, starts, ends, var_types

def parse_rsids(rsid_cells, num_variants):
    rsids = []
    for rsid_idx in range(RSID_COL, num_variants + 1):
        rsid = rsid_cells.loc[rsid_idx]
        if pd.notna(rsid):
            match = re.match(RSID_PATTERN, rsid)
            if match:
                rsid = match.group()
        else:
            rsid = ""
        rsids.append(rsid)

    return rsids

def parse_alleles(allele_cells, num_variants):
    allele_cells.reset_index(inplace=True, drop=True)
    #var_types =

    
    if allele_cells.iloc[:, 0].str.contains("NOTES").eq(False).all():
        hap_end_row = len(allele_cells.index)
    else:
        hap_end_row = allele_cells.index[allele_cells.iloc[:, 0] == "NOTES:"].values.item(0) - 1
   
    allele_cells = allele_cells[:hap_end_row]
    for col_idx in range(VARIANT_COL, num_variants + VARIANT_COL):
        ref_allele = allele_cells.loc[0, col_idx]
        allele_cells.iloc[:, col_idx] = allele_cells.iloc[:, col_idx].replace({np.nan: ref_allele})

        #var_type = "SNP"
        #allele_col = allele_cells.iloc[:, col_idx]
        #startswith_del_alleles = set(allele_col.loc[allele_col.str.startswith("del")])
        #if len(startswith_del_alleles) > 0:
            #for a in startswith_del_alleles:
                #if a != "delGene":
                    #if ref_allele == a:
                        #var_type = "INS"
                    #e
                    # lse: var_type = ("DEL")
                    #break
        #var_types.append(var_type)

    #####
    global alleles
    alleles = allele_cells.copy()

    allele_cells['alleles'] = allele_cells.values.tolist()
    allele_cells['alleles'] = allele_cells['alleles'].apply(lambda l: ','.join(l[2:]))
    haps = allele_cells.filter([0,'alleles'], axis=1)
    haps.columns = ["name", "alleles"]

    #####
    global hap_names
    hap_names = haps['name']

    return haps

def parse():
    files = filter(lambda f: f.endswith('.xlsx'), listdir(HAPLOTYPE_TABLE_DIR))
    for file in files:
        definition_file = path.join(HAPLOTYPE_TABLE_DIR, file)
        definition_table = pd.read_excel(definition_file, header=None)
        
        gene = parse_gene(definition_table.iloc[GENE_ROW][GENE_COL])

        print(gene)
        
        if gene == 'G6PD':
            global CHROM_COL, VARIANT_COL, RSID_COL
            CHROM_COL = CHROM_COL + 1
            VARIANT_COL = VARIANT_COL + 1
            RSID_COL = RSID_COL + 1

        chrom = parse_chrom(definition_table.iloc[CHROM_ROW][CHROM_COL])
        num_variants = definition_table.iloc[CHROM_ROW, VARIANT_COL:].count()
        chrom_hgvs_names, starts, ends, var_types = parse_variants(definition_table.iloc[VARIANT_ROW, VARIANT_COL:], num_variants)
        
        rsids = parse_rsids(definition_table.iloc[RSID_ROW, RSID_COL:num_variants+RSID_COL], num_variants)
        haps = parse_alleles(definition_table.iloc[HAP_ROW:, HAP_COL:num_variants + VARIANT_COL], num_variants)
        haps["gene"] = gene
        haps["chrom"] = chrom
        haps["num_variants"] = num_variants
        starts = ','.join(map(str, starts))
        haps["starts"] = starts
        ends = ','.join(map(str, ends))
        haps["ends"] = ends
        chrom_hgvs_names = ','.join(chrom_hgvs_names)
        haps["chrom_hgvs_names"] = chrom_hgvs_names
        rsids = ','.join(rsids)
        haps["rsids"] = rsids
        var_types = ','.join(var_types)
        haps["types"] = var_types
        # print(haps.loc[haps['name'] == "Hektoen"])
        haps.to_csv(path.join(OUT_DIR, gene+'.tsv'), sep='\t', index=False)
        
        
        ##########################For .json translation file #########################

        # variants
        global positions
        chrom_hgvs_names = chrom_hgvs_names.split(",")
        var_types = var_types.split(",")
        rsids = rsids.split(",")
        variants = []
        for pos, hgvs, vtype, rsid in zip(positions, chrom_hgvs_names, var_types, rsids):
            variants.append({
                "chromosome": chrom,
                "position": pos,
                "chrom_hgvs_name": hgvs,
                "type": vtype,
                "rsid": rsid
            })
        print("variants = " + str(len(variants)))
  
        # named_alleles
        global hap_names
        hap_names = hap_names.values.tolist()
        global alleles
        alleles = alleles.iloc[:, VARIANT_COL:].values.tolist()
        named_alleles = []

        for name, hapal in zip(hap_names, alleles):
            named_alleles.append({
                "name": name,
                "function": None,
                "alleles": dict(zip(positions, hapal))
            })

        json_definition = {
            "gene": gene,
            "chromosome": chrom,
            "variants": variants,
            "named_alleles": named_alleles
        }

        with open(path.join(JSON_DIR, gene + "_allele_definition_table.json"), 'w') as file:
            json.dump(json_definition, file, indent=4)
       
        if gene == "G6PD":
            CHROM_COL = CHROM_COL - 1
            VARIANT_COL = VARIANT_COL - 1
            RSID_COL = RSID_COL - 1

        
        positions = []

if __name__ == '__main__':
    parse()