#!/usr/bin/env python
from glob import glob
from os import path
import pandas as pd
import re
import numpy as np
from pprint import pprint
import json
import argparse

ALLELE_DEFINITION_DIR = path.join(path.dirname(__file__), '..', "data", "allele_definition_table")
DELIMITER_FILE = path.join(path.dirname(__file__), '..', "data", "allele_definition_table", "delimiter.tsv")
OUT_DIR = path.join(path.dirname(__file__), '..', "out") 

GENE_PATTERN = r'^GENE:\s*(\w+)$'
CHROM_PATTERN = r"chromosome\s(\w+)"
SNP_PATTERN = r'^g\.(\d+)([A-Z]>([A-Z]|[A-Z](\/[A-Z])*))*$'
INS_PATTERN = r'^g\.(\d+)_?(\d*)ins.+$'
DEL_PATTERN = r'^g\.(\d+)_?(\d*)del.+$'

# SV_DEL_PATTERN = r"^delGene$"
RSID_PATTERN = r'^rs\d+$'

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

def parse_gene(gene_cell):
	gene_cell = gene_cell.strip()

	match = re.match(GENE_PATTERN, gene_cell)
	if match:
		gene = match.group(1)
		return gene
	else:
		return

def parse_chrom(chrom_cell):
	chrom_name = None
	match = re.search(CHROM_PATTERN, chrom_cell)
	if match:
		chrom_name = "chr" + match.group(1)
		return chrom_name
	else:
		return

def parse_deletion(chrom_hgvs_name):
	match = re.match(DEL_PATTERN, chrom_hgvs_name)
	if match:
		start = match.group(1)
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

def parse_haps(gene, full_hap_cells, last_col):
	chrom_hgvs_names = [np.nan]
	starts = [np.nan]
	ends = [np.nan]
	var_types = [np.nan]

	for var_idx in range(1, last_col+1):
		chrom_hgvs_name = full_hap_cells.iloc[0, var_idx]
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
			print(f"Cannot determine variant type of {chrom_hgvs_name} of gene {gene}")

		chrom_hgvs_names.append(chrom_hgvs_name)
		starts.append(start)
		ends.append(end)
		
		ref_allele = full_hap_cells.loc[1, var_idx]
		full_hap_cells.iloc[1:, var_idx] = full_hap_cells.iloc[1:, var_idx].replace({np.nan: ref_allele})
		full_hap_cells.iloc[1:, var_idx] = full_hap_cells.iloc[1:, var_idx].apply(lambda allele: iupac[allele] if len(allele)==1 else allele)

	loci = pd.DataFrame([starts, ends, var_types])

	haps_by_loci = pd.concat([loci, full_hap_cells], ignore_index=True)
	haps_by_loci.loc[0] = haps_by_loci.loc[0].astype(float)
	haps_by_loci.loc[1] = haps_by_loci.loc[1].astype(float)

	_sorted = haps_by_loci.iloc[:, 1:].sort_values(by=0, axis=1)
	_combined_named_alleles = pd.concat([haps_by_loci[[0]], _sorted], axis=1)

	sorted_starts = _combined_named_alleles.iloc[0, 1:]
	sorted_ends = _combined_named_alleles.iloc[1, 1:]
	sorted_var_types = _combined_named_alleles.iloc[2, 1:]
	sorted_hgvs = _combined_named_alleles.iloc[3, 1:]
	named_alleles = _combined_named_alleles.iloc[4:, 0].values
	alleles = _combined_named_alleles.iloc[4:, 1:].values
	
	return sorted_starts, sorted_ends, sorted_var_types, sorted_hgvs, named_alleles, alleles

def parse_rsids(gene, rsid_cells, last_col):
	# print(gene)
	# print(rsid_cells)
	rsids = []
	for rsid_idx in range(0, last_col):
		rsid = rsid_cells.iloc[0, rsid_idx]
		# print(rsid_idx, rsid)
		if pd.notna(rsid):
			match = re.match(RSID_PATTERN, rsid)
			if match:
				rsid = match.group()
		else:
			rsid = ""
		rsids.append(rsid)
	return rsids

def parse():
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers(dest='subparser_name')

	parser_parse = subparsers.add_parser('parse')
	parser_parse.add_argument('--a', default=ALLELE_DEFINITION_DIR)
	parser_parse.add_argument('--d', default=DELIMITER_FILE)
	parser_parse.add_argument('--o', default=OUT_DIR)

	args = parser.parse_args()

	if args.subparser_name == 'parse':
		allele_dir = args.a
		delim_file = args.d
		out_dir = args.o

		allele_definitions = pd.DataFrame()

		delim = pd.read_csv(delim_file, sep='\t', index_col=0).to_dict('index')

		for fullpath in glob(f'{allele_dir}/*.xlsx'):
			definition_table = pd.read_excel(fullpath, header=None)

			keyword_gene = path.splitext(path.basename(fullpath))[0].split("_", 1)[0]

			last_col = delim[keyword_gene]['last_col']
			
			gene_row = delim[keyword_gene]['gene_row']
			gene_col = delim[keyword_gene]['gene_col']
			gene_cell = definition_table.iloc[gene_row][gene_col]
			gene = parse_gene(gene_cell)
		
			chrom_row = delim[keyword_gene]['chrom_row']
			chrom_col = delim[keyword_gene]['chrom_col']
			chrom_cell = definition_table.iloc[chrom_row][chrom_col]
			chrom = parse_chrom(chrom_cell)

			variant_row = delim[keyword_gene]['variant_row']
			variant_col = delim[keyword_gene]['variant_col']
			variant_cells = definition_table.iloc[variant_row, variant_col:last_col+1].to_frame().T

			hap_start_row = delim[keyword_gene]['hap_start_row']
			hap_end_row = delim[keyword_gene]['hap_end_row']
			hap_col = delim[keyword_gene]['hap_col']
			hap_cells = definition_table.iloc[hap_start_row:hap_end_row+1, hap_col:last_col+1]
			if gene == "G6PD":
				hap_cells = hap_cells.drop(1, axis=1)
				hap_cells = hap_cells.rename(columns=dict(zip(hap_cells.columns, range(0, last_col))))
				variant_cells = variant_cells.rename(columns=dict(zip(variant_cells.columns, range(1, last_col))))
				last_col = last_col - 1
	
			full_hap_cells = pd.concat([variant_cells, hap_cells], axis=0).reset_index(drop=True)

			starts, ends, var_types, hgvs, named_alleles, alleles = parse_haps(gene, full_hap_cells, last_col)

			rsid_row = delim[keyword_gene]['rsid_row']
			rsid_col = delim[keyword_gene]['rsid_col']
			rsid_cells = definition_table.iloc[rsid_row, rsid_col:last_col+1].to_frame().T
			rsid_cells = variant_cells.rename(columns=dict(zip(variant_cells.columns, range(0, last_col))))
			rsids = parse_rsids(gene, rsid_cells, last_col)

			# print(gene, starts.count(), ends.count(), var_types.count(), hgvs.count(), alleles[0].shape[0], len(rsids))
			assert starts.count() == ends.count() ==  var_types.count() == hgvs.count() == alleles[0].shape[0] == len(rsids)
			
			num_variants = starts.count()

		
			##################
			# Arrange tables #
			##################
			concat_alleles = [','.join(it) for it in alleles]
			round_starts = [round(it) for it in starts]
			round_ends = [round(it) for it in ends]

			allele_definitions = allele_definitions.append(pd.DataFrame(data={"name":named_alleles, 
				"gene":gene, 
				"chrom":chrom,
				"num_variants": num_variants,
				"starts": ','.join(map(str, round_starts)),
				"ends": ','.join(map(str, round_ends)),
				"chrom_hgvs_names": ','.join(map(str, hgvs)),
				"rsids": ','.join(map(str, rsids)),
				"var_types": ','.join(map(str, var_types)),
				"alleles": concat_alleles
				}), ignore_index=True)

		allele_definitions.to_csv(path.join(out_dir, 'allele_definitions.tsv'), sep='\t', index=False)

 
if __name__ == '__main__':
	parse()