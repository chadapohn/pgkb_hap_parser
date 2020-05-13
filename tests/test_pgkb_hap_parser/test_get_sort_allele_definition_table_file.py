#!/usr/bin/env python

import os
import unittest
from src.pgkb_hap_parser.get_sort_allele_definition_table_file import get_sort_allele_definition_table_file

def test_sort_allele_definition_file():
    allele_definition_path = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'allele_definition_table')
    expect = [
                allele_definition_path + '/CACNA1S_allele_definition_table.xlsx', 
                allele_definition_path + '/CFTR_allele_definition_table.xlsx', 
                allele_definition_path + '/CYP2B6_allele_definition_table.xlsx', 
                allele_definition_path + '/CYP2C9_allele_definition_table.xlsx', 
                allele_definition_path + '/CYP2C19_allele_definition_table.xlsx', 
                allele_definition_path + '/CYP2D6_allele_definition_table.xlsx', 
                allele_definition_path + '/CYP3A5_allele_definition_table.xlsx', 
                allele_definition_path + '/CYP4F2_allele_definition_table.xlsx', 
                allele_definition_path + '/DPYD_allele_definition_table.xlsx', 
                allele_definition_path + '/G6PD_allele_definition_table.xlsx', 
                allele_definition_path + '/IFNL3_allele_definition_table.xlsx', 
                allele_definition_path + '/NUDT15_allele_definition_table.xlsx', 
                allele_definition_path + '/RYR1_allele_definition_table.xlsx', 
                allele_definition_path + '/SLCO1B1_allele_definition_table.xlsx', 
                allele_definition_path + '/TPMT_allele_definition_table.xlsx', 
                allele_definition_path + '/UGT1A1_allele_definition_table.xlsx', 
                allele_definition_path + '/VKORC1_allele_definition_table.xlsx'
            ]
    actual = get_sort_allele_definition_table_file(allele_definition_path)
    assert expect == actual
