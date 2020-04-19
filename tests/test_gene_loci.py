from pandas.testing import assert_frame_equal

from pgkb_hap_parser.parser import parse_gene, parse_chrom, parse_variants

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




   
