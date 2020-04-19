#!/usr/bin/env python
from os import path, listdir
import pandas as pd
import re
import numpy as np
from pprint import pprint
import json

HAPLOTYPE_TABLE_DIR = path.join(path.dirname(__file__), '..', "data", "allele_definition_table")
OUT_DIR = path.join(path.dirname(__file__), '..', "out") 

######
GENE_ROW = 0
GENE_COL = 0
GENE_PATTERN = r'^GENE:\s*(\w+)$'

CHROM_ROW = 3
CHROM_COL = 0
CHROM_PATTERN = r"chromosome\s(\w+)"

VARIANT_ROW = 3
VARIANT_COL = 1
SNP_PATTERN = r'^g\.(\d+)[A-Z]>([A-Z]|[ATCG](\/[A-Z])*)$'
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
	match = re.search(CHROM_PATTERN, chrom_cell)

	if match:
		chrom_name = "chr" + match.group(1)
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
		end = match.group(1)

		return start, end
	else:
		return

def parse_variants(variant_cells, num_variants, VARIANT_COL=VARIANT_COL):
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
			print(10*"=")
			print(chrom_hgvs_name)
			start = "-"
			end = "-"

		starts.append(start)
		ends.append(end)

	#_sorted_variants = sorted(zip(starts, ends, chrom_hgvs_names, var_types), key=lambda x: int(x[0]))
	# starts, ends, chrom_hgvs_names, var_types = zip(*_sorted_variants)
	return chrom_hgvs_names, starts, ends, var_types

def parse_rsids(rsid_cells, num_variants):
    rsids = []
    for rsid_idx in range(RSID_COL, num_variants):
        rsid = rsid_cells.loc[rsid_idx]
        if pd.notna(rsid):
            match = re.match(RSID_PATTERN, rsid)
            if match:
                rsid = match.group()
        else:
            rsid = ""
        rsids.append(rsid)
    return rsids

def parse_alleles(gene, allele_cells, num_variants):
	if gene == 'G6PD':
		global VARIANT_COL
		VARIANT_COL = 2

	allele_cells.reset_index(inplace=True, drop=True)
	
	if allele_cells.iloc[:, 0].str.contains("NOTES").eq(False).all():
		hap_end_row = len(allele_cells.index)
	else:
		hap_end_row = allele_cells.index[allele_cells.iloc[:, 0] == "NOTES:"].values.item(0) - 1
   
	allele_cells = allele_cells[:hap_end_row]
	for col_idx in range(VARIANT_COL, num_variants + VARIANT_COL):
		ref_allele = allele_cells.loc[0, col_idx]
		allele_cells.iloc[:, col_idx] = allele_cells.iloc[:, col_idx].replace({np.nan: ref_allele})
  
	#####
	global alleles
	alleles = allele_cells.copy()
	
	allele_cells['alleles'] = allele_cells.values.tolist()	
	allele_cells['alleles'] = allele_cells['alleles'].apply(lambda l: ','.join(map(str, l[VARIANT_COL:])))
	haps = allele_cells.filter([0,'alleles'], axis=1)
	haps.columns = ["name", "alleles"]

	#####
	global hap_names
	hap_names = haps['name']
	return haps, hap_end_row

def parse():
	table = pd.DataFrame()
	files = filter(lambda f: f.endswith('.xlsx'), listdir(HAPLOTYPE_TABLE_DIR))
	for file in files:
		definition_file = path.join(HAPLOTYPE_TABLE_DIR, file)
		definition_table = pd.read_excel(definition_file, header=None)
		definition_table = definition_table.sort_values(axis=1, by=3, na_position='first')
		definition_table.columns = range(definition_table.shape[1])
		
		gene = parse_gene(definition_table.iloc[GENE_ROW][GENE_COL])
		
		if gene == 'G6PD':
			global CHROM_COL, VARIANT_COL, RSID_COL
			CHROM_COL = 1
			VARIANT_COL = 2
			RSID_COL = 2

		gchrom = parse_chrom(definition_table.iloc[CHROM_ROW][CHROM_COL])
		num_variants = definition_table.iloc[CHROM_ROW, VARIANT_COL:].count()
		chrom_hgvs_names, starts, ends, var_types = parse_variants(definition_table.iloc[VARIANT_ROW, VARIANT_COL:], num_variants, VARIANT_COL)
		if gene == 'G6PD':
			rsids = parse_rsids(definition_table.iloc[RSID_ROW, RSID_COL:num_variants+RSID_COL], num_variants + 2)
		else:
			rsids = parse_rsids(definition_table.iloc[RSID_ROW, RSID_COL:num_variants+RSID_COL], num_variants + 1)
		haps, hap_end_row = parse_alleles(gene, definition_table.iloc[HAP_ROW:, HAP_COL:num_variants + VARIANT_COL], num_variants)

		# if gene != "G6PD":
		# 	print(gene)
		# 	cells = definition_table.iloc[HAP_ROW:, HAP_COL:num_variants + VARIANT_COL]
		
		# 	# cells = cells.drop(1, axis=1) G6PD only
		# 	cells = cells.set_index(0)
		# 	cells.columns = chrom_hgvs_names
		# 	for idx in range(1, len(cells.index)):
		# 		cell = cells.iloc[idx]
		# 		res = {}
		# 		res[cell.name] = []
		# 		for allele, hgvs in zip(cell.values, cell.index):
		# 			if allele is not np.nan:
		# 				res[cell.name].append(hgvs)
		# 		print(cell.name, ', '.join(res[cell.name]).encode('utf-8'))

		print(gene)
		cells = definition_table.iloc[HAP_ROW + 1:, HAP_COL:num_variants + VARIANT_COL]
		if gene == "G6PD":
			cells = cells.drop(1, axis=1)
		cells = cells.set_index(0)
		cells.columns = chrom_hgvs_names
		cells = cells.T
		for idx in range(0, len(cells.index)):
			cell = cells.iloc[idx]
			res = {}
			res[cell.name] = []
			for allele, hgvs in zip(cell.values, cell.index):
				if allele is not np.nan:
					res[cell.name].append(hgvs)
			print(cell.name, ', '.join(res[cell.name]).encode('utf-8'))
	
		##################
		# Arrange tables #
		##################
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

	table.to_csv(path.join(OUT_DIR, 'allele.definitions.tsv'), sep='\t', index=False)

 
if __name__ == '__main__':
	parse()