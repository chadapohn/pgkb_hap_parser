from pandas.testing import assert_frame_equal

from pgkb_hap_parser.parser import parse_gene, parse_chrom

def test_parse_gene():
    assert parse_gene('GENE:CACNA1S') == 'CACNA1S'
    assert parse_gene('GENE:CFTR') == 'CFTR'
    assert parse_gene('GENE: CYP2B6') == 'CYP2B6'
    assert parse_gene('GENE:CYP2C9') == 'CYP2C9'
    assert parse_gene('GENE: CYP2C19 ') == 'CYP2C19'
    assert parse_gene('GENE:CYP3A5') == 'CYP3A5'
    assert parse_gene('GENE:CYP4F2') == 'CYP4F2'
    assert parse_gene('GENE:DPYD') == 'DPYD'
    assert parse_gene('GENE:G6PD') == 'G6PD'
    assert parse_gene('GENE:IFNL3') == 'IFNL3'
    assert parse_gene('GENE: NUDT15') == 'NUDT15'
    assert parse_gene('GENE:RYR1') == 'RYR1'
    assert parse_gene('GENE: SLCO1B1') == 'SLCO1B1'
    assert parse_gene('GENE: TPMT') == 'TPMT'
    assert parse_gene('GENE:VKORC1') == 'VKORC1'

def test_parse_chrom():
    assert parse_chrom('Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)') == 'chr1'
    assert parse_chrom('Position at NC_000013.11 (Homo sapiens chromosome 6, GRCh38.p7)') == 'chr6'
    assert parse_chrom('Position at NC_000007.14 (Homo sapiens chromosome 7, GRCh38.p2)') == 'chr7'
    assert parse_chrom('Position at NC_000010.11 (Homo sapiens chromosome 10, GRCh38.p2)') == 'chr10'
    assert parse_chrom('Position at NC_000012.12 (Homo sapiens chromosome 12, GRCh38.p2)') == 'chr12'
    assert parse_chrom('Position at NC_000019.10 (Homo sapiens chromosome 16, GRCh38.p2)') == 'chr16'
    assert parse_chrom('Position at NC_000019.10 (Homo sapiens chromosome 19, GRCh38.p2)') == 'chr19'
    assert parse_chrom('Position at NC_000023.11 (Homo sapiens chromosome X, GRCh38.p2)') == 'chrX'

# def test_parse_variants():
#     CHROM_ROW = 3
#     CHROM_COL, VARIANT_COL = 0, 1
#     RSID_ROW = 5
#     RSID_COL = 1

#     CACNA1S = os.getcwd() + '/data/xlsx/CACNA1S_allele_definition_table.xlsx'
#     df = pd.read_excel(CACNA1S, header=None)
#     num_variants = df.iloc[CHROM_ROW, VARIANT_COL:].count()
#     actual_chrom_hgvs_names, actual_starts, actual_ends, actual_var_types = parse_variants(df.iloc[CHROM_ROW, VARIANT_COL:], num_variants)
#     expect_chrom_hgvs_names = ['g.201091993G>A', 'g.201060815C>T']
#     expect_starts = ['201091992', '201060814']
#     expect_ends = ['201091993', '201060815']
#     expect_var_types = ['SNP', 'SNP']
#     assert actual_chrom_hgvs_names == expect_chrom_hgvs_names
#     assert actual_starts == expect_starts
#     assert actual_ends == expect_ends
#     assert actual_var_types == expect_var_types

#     VARIANT_COL = 13

#     CYP2C9 = os.getcwd() + '/data/xlsx/CYP2C9_allele_definition_table.xlsx'
#     df = pd.read_excel(CYP2C9, header=None)
#     actual_chrom_hgvs_names, actual_starts, actual_ends, actual_var_types = parse_variants(df.iloc[CHROM_ROW, VARIANT_COL:VARIANT_COL + 1], 1, VARIANT_COL)
#     expect_chrom_hgvs_names = ['g.94942213_94942222delAGAAATGGAA;g.94942216A>G']
#     expect_starts = ['94942212']
#     expect_ends = ['94942222']
#     expect_var_types = ['DEL']
#     assert actual_chrom_hgvs_names == expect_chrom_hgvs_names
#     assert actual_starts == expect_starts
#     assert actual_ends == expect_ends
#     assert actual_var_types == expect_var_types

#     VARIANT_COL = 32

#     CYP2C9 = os.getcwd() + '/data/xlsx/CYP2C9_allele_definition_table.xlsx'
#     df = pd.read_excel(CYP2C9, header=None)
#     actual_chrom_hgvs_names, actual_starts, actual_ends, actual_var_types = parse_variants(df.iloc[CHROM_ROW, VARIANT_COL:VARIANT_COL + 1], 1, VARIANT_COL)
#     expect_chrom_hgvs_names = ['g.94949283delA']
#     expect_starts = ['94949282']
#     expect_ends = ['94949283']
#     expect_var_types = ['DEL']
#     assert actual_chrom_hgvs_names == expect_chrom_hgvs_names
#     assert actual_starts == expect_starts
#     assert actual_ends == expect_ends
#     assert actual_var_types == expect_var_types

#     VARIANT_COL = 5

#     NUDT15 = os.getcwd() + '/data/xlsx/NUDT15_allele_definition_table.xlsx'
#     df = pd.read_excel(NUDT15, header=None)
#     actual_chrom_hgvs_names, actual_starts, actual_ends, actual_var_types = parse_variants(df.iloc[CHROM_ROW, VARIANT_COL:VARIANT_COL + 1], 1, VARIANT_COL)
#     expect_chrom_hgvs_names = ['g.48037801_48037802insGAGTCG']
#     expect_starts = ['48037801']
#     expect_ends = ['48037802']
#     expect_var_types = ['INS']
#     assert actual_chrom_hgvs_names == expect_chrom_hgvs_names
#     assert actual_starts == expect_starts
#     assert actual_ends == expect_ends
#     assert actual_var_types == expect_var_types

#     VARIANT_COL = 13

#     NUDT15 = os.getcwd() + '/data/xlsx/NUDT15_allele_definition_table.xlsx'
#     df = pd.read_excel(NUDT15, header=None)
#     actual_chrom_hgvs_names, actual_starts, actual_ends, actual_var_types = parse_variants(df.iloc[CHROM_ROW, VARIANT_COL:VARIANT_COL + 1], 1, VARIANT_COL)
#     expect_chrom_hgvs_names = ['g.48041103_48041104insG']
#     expect_starts = ['48041103']
#     expect_ends = ['48041104']
#     expect_var_types = ['INS']
#     assert actual_chrom_hgvs_names == expect_chrom_hgvs_names
#     assert actual_starts == expect_starts
#     assert actual_ends == expect_ends
#     assert actual_var_types == expect_var_types

# def test_parse_rsid():
#     CHROM_ROW = 3
#     CHROM_COL, VARIANT_COL = 0, 1
#     RSID_ROW = 5
#     RSID_COL = 1

#     CACNA1S = os.getcwd() + '/data/xlsx/CACNA1S_allele_definition_table.xlsx'
#     df = pd.read_excel(CACNA1S, header=None)
#     num_variants = df.iloc[CHROM_ROW, VARIANT_COL:].count()
#     actual = parse_rsids(df.iloc[RSID_ROW, RSID_COL:RSID_COL + num_variants], num_variants)
#     expect = ['rs772226819', 'rs1800559']
#     assert actual == expect

#     CYP4F2 = os.getcwd() + '/data/xlsx/CYP4F2_allele_definition_table.xlsx'
#     df = pd.read_excel(CYP4F2, header=None)
#     num_variants = df.iloc[CHROM_ROW, VARIANT_COL:].count()
#     actual = parse_rsids(df.iloc[RSID_ROW, RSID_COL:RSID_COL + num_variants], num_variants)
#     expect = ['rs3093105', 'rs2108622']
#     assert actual == expect

#     UGT1A1 = os.getcwd() + '/data/xlsx/UGT1A1_allele_definition_table.xlsx'
#     df = pd.read_excel(UGT1A1, header=None)
#     num_variants = df.iloc[CHROM_ROW, VARIANT_COL:].count()
#     actual = parse_rsids(df.iloc[RSID_ROW, RSID_COL:RSID_COL + num_variants], num_variants)
#     expect = ['rs4124874', 'rs887829', '', 'rs4148323', 'rs35350960']
#     assert actual == expect