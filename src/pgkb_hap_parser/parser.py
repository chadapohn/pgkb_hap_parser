#!/usr/bin/env python

import time
import argparse
from src.pgkb_hap_parser.get_sort_allele_definition_table_file import get_sort_allele_definition_table_file
import pandas as pd
from src.pgkb_hap_parser.get_gene import get_gene
from src.pgkb_hap_parser.get_chromosome import get_chromosome
from src.pgkb_hap_parser.get_hgvs_start_end_variant_type import get_hgvs_start_end_variant_type
from src.pgkb_hap_parser.get_rsid import get_rsid
from src.pgkb_hap_parser.get_haplotype_name_variant_extract_iupac import get_haplotype_name_variant_extract_iupac
from src.pgkb_hap_parser.get_name_relation_to_hgvs import get_name_relation_to_hgvs
from src.pgkb_hap_parser.get_hgvs_relation_to_name import get_hgvs_relation_to_name

def parse():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subparser_name")

    parser_parse = subparsers.add_parser("parse")
    parser_parse.add_argument("-d", "--delimiter", help="Path to the delimiter file", required=True)
    parser_parse.add_argument("-a", "--allele_definition_table", help="Path to the allele definition table directory", required=True)
    parser_parse.add_argument("-o", "--output", help="Path to the output directory", required=True)
    
    args = parser.parse_args()

    if args.subparser_name == "parse":
        start_time = time.time()

        delimiter_path = args.delimiter
        allele_definition_table_path = args.allele_definition_table
        output_path = args.output

        delimiter_df = pd.read_csv(delimiter_path, index_col="gene", sep="\t")

        sort_allele_definition_table_file = get_sort_allele_definition_table_file(allele_definition_table_path)
        
        allele_definition_df = pd.DataFrame(columns=["name", "gene", "chromosome", "hgvs", "start", "end", "rsid", "variant_type", "variant"])
        allele_definition_name_relation_to_hgvs_df = pd.DataFrame(columns=["gene", "name", "hgvs"])
        allele_definition_hgvs_relation_to_name_df = pd.DataFrame(columns=["gene", "hgvs", "name"])

        for allele_definition_table in sort_allele_definition_table_file:
            allele_definition_table_df = pd.read_excel(allele_definition_table, header=None)
            
            gene_row = 0
            gene_col = 0
            gene_cell = allele_definition_table_df.iloc[gene_row, gene_col]
            gene = get_gene(gene_cell)

            try:
                chromosome_row = delimiter_df.loc[gene, "chromosome_row"]
                chromosome_col = delimiter_df.loc[gene, "chromosome_col"]
                chromosome_cell = allele_definition_table_df.iloc[chromosome_row, chromosome_col]
                chromosome = get_chromosome(chromosome_cell)
            except:
                print(f"not have information about gene \"{gene}\" in delimiter.tsv, please add information before and try it again.")
                exit()

            allele_definition_table_df.sort_values(axis=1, by=int(delimiter_df.loc[gene, "hgvs_row"]), inplace=True)
            allele_definition_table_df.columns = range(allele_definition_table_df.shape[1])

            hgvs_row = delimiter_df.loc[gene, "hgvs_row"]
            hgvs_start_col = delimiter_df.loc[gene, "hgvs_start_col"]
            hgvs_end_col = delimiter_df.loc[gene, "hgvs_end_col"]
            hgvs_cell = allele_definition_table_df.iloc[hgvs_row, hgvs_start_col:hgvs_end_col]
            hgvs, start, end, variant_type = get_hgvs_start_end_variant_type(hgvs_cell)

            rsid_row = delimiter_df.loc[gene, "rsid_row"]
            rsid_start_col = delimiter_df.loc[gene, "rsid_start_col"]
            rsid_end_col = delimiter_df.loc[gene, "rsid_end_col"]
            rsid_cell = allele_definition_table_df.iloc[rsid_row, rsid_start_col:rsid_end_col]
            rsid = get_rsid(rsid_cell)

            haplotype_start_row = delimiter_df.loc[gene, "haplotype_start_row"]
            haplotype_end_row = delimiter_df.loc[gene, "haplotype_end_row"]
            haplotype_start_col = delimiter_df.loc[gene, "haplotype_start_col"]
            haplotype_end_col = delimiter_df.loc[gene, "haplotype_end_col"]
            haplotype_cell = allele_definition_table_df.iloc[haplotype_start_row:haplotype_end_row, haplotype_start_col:haplotype_end_col]
            haplotype_extract_iupac, name_extract_iupac, variant_extract_iupac = get_haplotype_name_variant_extract_iupac(haplotype_cell)
            
            assert len(hgvs) == len(start) == len(end) == len(variant_type) == len(rsid) == len(variant_extract_iupac[0])
            assert len(haplotype_extract_iupac) == len(name_extract_iupac) == len(variant_extract_iupac)
            
            allele_definition_extract_iupac_df = pd.DataFrame(list(zip([gene], [chromosome], [hgvs], [start], [end], [rsid], [variant_type])), columns=["gene", "chromosome", "hgvs", "start", "end", "rsid", "variant_type"])
            allele_definition_extract_iupac_df = allele_definition_extract_iupac_df.append([allele_definition_extract_iupac_df] * (len(haplotype_extract_iupac)- 1), ignore_index=True)
            allele_definition_extract_iupac_df.insert(0, "name", name_extract_iupac)
            allele_definition_extract_iupac_df.insert(allele_definition_extract_iupac_df.shape[1], "variant", variant_extract_iupac)
            allele_definition_df = pd.concat([allele_definition_df, allele_definition_extract_iupac_df], ignore_index=True)

            name_relation_to_hgvs = get_name_relation_to_hgvs(hgvs, haplotype_cell)
            name_relation_to_hgvs_df = pd.DataFrame([gene], columns=["gene"])
            name_relation_to_hgvs_df = name_relation_to_hgvs_df.append([name_relation_to_hgvs_df] * (len(name_relation_to_hgvs.keys()) - 1))
            name_relation_to_hgvs_df.insert(name_relation_to_hgvs_df.shape[1], "name", name_relation_to_hgvs.keys())
            name_relation_to_hgvs_df.insert(name_relation_to_hgvs_df.shape[1], "hgvs", list(name_relation_to_hgvs.values()))
            allele_definition_name_relation_to_hgvs_df = pd.concat([allele_definition_name_relation_to_hgvs_df, name_relation_to_hgvs_df], ignore_index=True)

            hgvs_relation_to_name = get_hgvs_relation_to_name(hgvs, haplotype_cell)
            hgvs_relation_to_name_df = pd.DataFrame([gene], columns=["gene"])
            hgvs_relation_to_name_df = hgvs_relation_to_name_df.append([hgvs_relation_to_name_df] * (len(hgvs_relation_to_name.keys()) - 1))
            hgvs_relation_to_name_df.insert(hgvs_relation_to_name_df.shape[1], "hgvs", hgvs_relation_to_name.keys())
            hgvs_relation_to_name_df.insert(hgvs_relation_to_name_df.shape[1], "name", list(hgvs_relation_to_name.values()))
            allele_definition_hgvs_relation_to_name_df = pd.concat([allele_definition_hgvs_relation_to_name_df, hgvs_relation_to_name_df], ignore_index=True)
        
        allele_definition_df.to_csv(output_path + "/allele_definitions.tsv", index=False, sep="\t")
        allele_definition_name_relation_to_hgvs_df.to_csv(output_path + "/allele_definitions_name_relation_to_hgvs.tsv", index=False, sep="\t")
        allele_definition_hgvs_relation_to_name_df.to_csv(output_path + "/allele_definitions_hgvs_relation_to_name.tsv", index=False, sep="\t")

        print(f"run pgkb_hap_parser package successfully in {time.time() - start_time:.2f} seconds")

if __name__ == "__main__":
    parse()
