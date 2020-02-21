#!/usr/bin/env python
from os import path, listdir
import pandas as pd
import re
import numpy as np
from pprint import pprint
import json

HAPLOTYPE_TABLE_DIR = path.join(path.dirname(__file__), '..', "data", "xlsx")
OUT_DIR = path.join(path.dirname(__file__), '..', "out") 

######
GENE_ROW = 0
GENE_COL = 0
GENE_PATTERN = "^GENE:\s*(\w+)$"

CHROM_ROW = 3
CHROM_COL = 0
CHROM_PATTERN = r"^.*N\w_\s*(\d+)\.\d+"

VARIANT_ROW = 3
VARIANT_COL = 1
SNP_PATTERN = r'^g\.(\d+)[ATCG]>[ATCG]$'
INS_PATTERN = r'^g\.(\d+)_?(\d*)ins.+$'
DEL_PATTERN = r'^g\.(\d+)_?(\d*)del.+$'
# REPEAT_PATTERN = r"^[cgp]\.(\d+)"
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

def parse_deletion(chrom_hgvs_name):
	match = re.match(DEL_PATTERN, chrom_hgvs_name)
	if match:
		start = match.group(1)
		start = int(start) - 1
		start = str(start)
		if match.group(1) and not match.group(2):
			end = match.group(1)
		elif match.group(1) and match.group(2):
			end = match.group(2)
	
		return start, end
	else:
		return 

def parse_insertion(chrom_hgvs_name):
	match = re.match(INS_PATTERN, chrom_hgvs_name)
	if match:
		start = match.group(1)
		end = match.group(2)
	
		return start, end
	else:
		return

def parse_snp(chrom_hgvs_name):
	match = re.match(SNP_PATTERN, chrom_hgvs_name)
	if match:
		start = match.group(1)
		start = int(start) - 1
		start = str(start)
		end = match.group(1)

		return start, end
	else:
		return

def parse_variants(variant_cells, num_variants, VARIANT_COL=VARIANT_COL):
	print(CHROM_COL)
	print(VARIANT_COL)
	print(RSID_COL)
	chrom_hgvs_names = []
	starts = []
	ends = []
	var_types = []

	for var_idx in range(VARIANT_COL, num_variants + VARIANT_COL):
		chrom_hgvs_name = variant_cells.loc[var_idx]
		chrom_hgvs_names.append(chrom_hgvs_name)

		if parse_deletion(chrom_hgvs_name):
			start, end = parse_deletion(chrom_hgvs_name)
			var_types.append('DEL')

		elif parse_insertion(chrom_hgvs_name):
			start, end = parse_insertion(chrom_hgvs_name)
			var_types.append('INS')

		elif parse_snp(chrom_hgvs_name):
			start, end = parse_snp(chrom_hgvs_name)
			var_types.append('SNP')

		else:
			start = "-"
			end = "-"

		starts.append(start)
		ends.append(end)

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
	allele_cells['alleles'] = allele_cells['alleles'].apply(lambda l: ','.join(l[VARIANT_COL:]))
	haps = allele_cells.filter([0,'alleles'], axis=1)
	haps.columns = ["name", "alleles"]

	#####
	global hap_names
	hap_names = haps['name']

	return haps

def parse():
	table = pd.DataFrame()
	files = filter(lambda f: f.endswith('.xlsx'), listdir(HAPLOTYPE_TABLE_DIR))
	for file in files:
		definition_file = path.join(HAPLOTYPE_TABLE_DIR, file)
		definition_table = pd.read_excel(definition_file, header=None)
		
		gene = parse_gene(definition_table.iloc[GENE_ROW][GENE_COL])

		print(gene)
		
		if gene == 'G6PD':
			global CHROM_COL, VARIANT_COL, RSID_COL
			CHROM_COL = 1
			VARIANT_COL = 2
			RSID_COL = 2

		chrom = parse_chrom(definition_table.iloc[CHROM_ROW][CHROM_COL])
		num_variants = definition_table.iloc[CHROM_ROW, VARIANT_COL:].count()
		chrom_hgvs_names, starts, ends, var_types = parse_variants(definition_table.iloc[VARIANT_ROW, VARIANT_COL:], num_variants, VARIANT_COL)
		
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
		table = table.append(haps, ignore_index=True)

		
		CHROM_COL = 0
		VARIANT_COL = 1
		RSID_COL = 1
	table.to_csv(path.join(OUT_DIR, 'allele_definition.tsv'), sep='\t', index=False)

 
if __name__ == '__main__':
	parse()