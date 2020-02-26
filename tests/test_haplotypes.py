from pgkb_hap_parser.parser import *

import os
from pandas._testing import assert_frame_equal

def test_parse_gene():
    CHROM_ROW = 3
    CHROM_COL, VARIANT_COL = 0, 1

    CACNA1S = os.getcwd() + '/data/xlsx/CACNA1S_allele_definition_table.xlsx'
    df = pd.read_excel(CACNA1S, header=None)
    actual = parse_gene(df.iloc[GENE_ROW][GENE_COL])
    expect = 'CACNA1S'
    assert actual == expect

    UGT1A1 = os.getcwd() + '/data/xlsx/UGT1A1_allele_definition_table.xlsx'
    df = pd.read_excel(UGT1A1, header=None)
    actual = parse_gene(df.iloc[GENE_ROW][GENE_COL])
    expect = 'UGT1A1'
    assert actual == expect

    CHROM_ROW = 3
    CHROM_COL, VARIANT_COL = 1, 2

    G6PD = os.getcwd() + '/data/xlsx/G6PD_allele_definition_table.xlsx'
    df = pd.read_excel(G6PD, header=None)
    actual = parse_gene(df.iloc[GENE_ROW][GENE_COL])
    expect = 'G6PD'
    assert actual == expect

def test_parse_chrome():
    CHROM_ROW = 3
    CHROM_COL, VARIANT_COL = 0, 1

    CACNA1S = os.getcwd() + '/data/xlsx/CACNA1S_allele_definition_table.xlsx'
    df = pd.read_excel(CACNA1S, header=None)
    actual = parse_chrom(df.iloc[CHROM_ROW][CHROM_COL])
    expect = 'chr1'
    assert actual == expect

    UGT1A1 = os.getcwd() + '/data/xlsx/UGT1A1_allele_definition_table.xlsx'
    df = pd.read_excel(UGT1A1, header=None)
    actual = parse_chrom(df.iloc[CHROM_ROW][CHROM_COL])
    expect = 'chr2'
    assert actual == expect

    CHROM_ROW = 3
    CHROM_COL, VARIANT_COL = 1, 2

    G6PD = os.getcwd() + '/data/xlsx/G6PD_allele_definition_table.xlsx'
    df = pd.read_excel(G6PD, header=None)
    actual = parse_chrom(df.iloc[CHROM_ROW][CHROM_COL])
    expect = 'chrX'
    assert actual == expect

def test_parse_variants():
    CHROM_ROW = 3
    CHROM_COL, VARIANT_COL = 0, 1
    RSID_ROW = 5
    RSID_COL = 1

    CACNA1S = os.getcwd() + '/data/xlsx/CACNA1S_allele_definition_table.xlsx'
    df = pd.read_excel(CACNA1S, header=None)
    num_variants = df.iloc[CHROM_ROW, VARIANT_COL:].count()
    actual_chrom_hgvs_names, actual_starts, actual_ends, actual_var_types = parse_variants(df.iloc[CHROM_ROW, VARIANT_COL:], num_variants)
    expect_chrom_hgvs_names = ['g.201091993G>A', 'g.201060815C>T']
    expect_starts = ['201091992', '201060814']
    expect_ends = ['201091993', '201060815']
    expect_var_types = ['SNP', 'SNP']
    assert actual_chrom_hgvs_names == expect_chrom_hgvs_names
    assert actual_starts == expect_starts
    assert actual_ends == expect_ends
    assert actual_var_types == expect_var_types

    VARIANT_COL = 13

    CYP2C9 = os.getcwd() + '/data/xlsx/CYP2C9_allele_definition_table.xlsx'
    df = pd.read_excel(CYP2C9, header=None)
    actual_chrom_hgvs_names, actual_starts, actual_ends, actual_var_types = parse_variants(df.iloc[CHROM_ROW, VARIANT_COL:VARIANT_COL + 1], 1, VARIANT_COL)
    expect_chrom_hgvs_names = ['g.94942213_94942222delAGAAATGGAA;g.94942216A>G']
    expect_starts = ['94942212']
    expect_ends = ['94942222']
    expect_var_types = ['DEL']
    assert actual_chrom_hgvs_names == expect_chrom_hgvs_names
    assert actual_starts == expect_starts
    assert actual_ends == expect_ends
    assert actual_var_types == expect_var_types

    VARIANT_COL = 32

    CYP2C9 = os.getcwd() + '/data/xlsx/CYP2C9_allele_definition_table.xlsx'
    df = pd.read_excel(CYP2C9, header=None)
    actual_chrom_hgvs_names, actual_starts, actual_ends, actual_var_types = parse_variants(df.iloc[CHROM_ROW, VARIANT_COL:VARIANT_COL + 1], 1, VARIANT_COL)
    expect_chrom_hgvs_names = ['g.94949283delA']
    expect_starts = ['94949282']
    expect_ends = ['94949283']
    expect_var_types = ['DEL']
    assert actual_chrom_hgvs_names == expect_chrom_hgvs_names
    assert actual_starts == expect_starts
    assert actual_ends == expect_ends
    assert actual_var_types == expect_var_types

    VARIANT_COL = 5

    NUDT15 = os.getcwd() + '/data/xlsx/NUDT15_allele_definition_table.xlsx'
    df = pd.read_excel(NUDT15, header=None)
    actual_chrom_hgvs_names, actual_starts, actual_ends, actual_var_types = parse_variants(df.iloc[CHROM_ROW, VARIANT_COL:VARIANT_COL + 1], 1, VARIANT_COL)
    expect_chrom_hgvs_names = ['g.48037801_48037802insGAGTCG']
    expect_starts = ['48037801']
    expect_ends = ['48037802']
    expect_var_types = ['INS']
    assert actual_chrom_hgvs_names == expect_chrom_hgvs_names
    assert actual_starts == expect_starts
    assert actual_ends == expect_ends
    assert actual_var_types == expect_var_types

    VARIANT_COL = 13

    NUDT15 = os.getcwd() + '/data/xlsx/NUDT15_allele_definition_table.xlsx'
    df = pd.read_excel(NUDT15, header=None)
    actual_chrom_hgvs_names, actual_starts, actual_ends, actual_var_types = parse_variants(df.iloc[CHROM_ROW, VARIANT_COL:VARIANT_COL + 1], 1, VARIANT_COL)
    expect_chrom_hgvs_names = ['g.48041103_48041104insG']
    expect_starts = ['48041103']
    expect_ends = ['48041104']
    expect_var_types = ['INS']
    assert actual_chrom_hgvs_names == expect_chrom_hgvs_names
    assert actual_starts == expect_starts
    assert actual_ends == expect_ends
    assert actual_var_types == expect_var_types

def test_parse_rsid():
    CHROM_ROW = 3
    CHROM_COL, VARIANT_COL = 0, 1
    RSID_ROW = 5
    RSID_COL = 1

    CACNA1S = os.getcwd() + '/data/xlsx/CACNA1S_allele_definition_table.xlsx'
    df = pd.read_excel(CACNA1S, header=None)
    num_variants = df.iloc[CHROM_ROW, VARIANT_COL:].count()
    actual = parse_rsids(df.iloc[RSID_ROW, RSID_COL:RSID_COL + num_variants], num_variants)
    expect = ['rs772226819', 'rs1800559']
    assert actual == expect

    CYP4F2 = os.getcwd() + '/data/xlsx/CYP4F2_allele_definition_table.xlsx'
    df = pd.read_excel(CYP4F2, header=None)
    num_variants = df.iloc[CHROM_ROW, VARIANT_COL:].count()
    actual = parse_rsids(df.iloc[RSID_ROW, RSID_COL:RSID_COL + num_variants], num_variants)
    expect = ['rs3093105', 'rs2108622']
    assert actual == expect

    UGT1A1 = os.getcwd() + '/data/xlsx/UGT1A1_allele_definition_table.xlsx'
    df = pd.read_excel(UGT1A1, header=None)
    num_variants = df.iloc[CHROM_ROW, VARIANT_COL:].count()
    actual = parse_rsids(df.iloc[RSID_ROW, RSID_COL:RSID_COL + num_variants], num_variants)
    expect = ['rs4124874', 'rs887829', '', 'rs4148323', 'rs35350960']
    assert actual == expect

def test_parse_alleles():
    CHROM_ROW = 3
    CHROM_COL, VARIANT_COL = 0, 1
    HAP_ROW = 7
    HAP_COL = 0

    CACNA1S = os.getcwd() + '/data/xlsx/CACNA1S_allele_definition_table.xlsx'
    df = pd.read_excel(CACNA1S, header=None)
    num_variants = df.iloc[CHROM_ROW, VARIANT_COL:].count()
    actual = parse_alleles(df.iloc[HAP_ROW:, HAP_COL:VARIANT_COL + num_variants], num_variants)
    name = ['Reference', 'c.520C>T', 'c.3257G>A']
    alleles = ['G,C', 'A,C', 'G,T']
    expect = df = pd.DataFrame(list(zip(name, alleles)), columns =['name', 'alleles'])
    assert_frame_equal(actual, expect, check_dtype=False)

    # IFNL3 = os.getcwd() + '/data/xlsx/IFNL3_allele_definition_table.xlsx'
    # df = pd.read_excel(IFNL3, header=None)
    # num_variants = df.iloc[CHROM_ROW, VARIANT_COL:].count()
    # actual = parse_alleles(df.iloc[HAP_ROW:, HAP_COL:VARIANT_COL + num_variants], num_variants)
    # name = ['rs12979860 reference (C)', 'rs12979860 variant (T)']
    # alleles = ['C', 'T']
    # expect = df = pd.DataFrame(list(zip(name, alleles)), columns =['name', 'alleles'])
    # assert_frame_equal(actual, expect, check_dtype=False)
