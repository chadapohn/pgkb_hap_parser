from pharmgkb-haplotype-table-parser.parser import parse_gene

def test_parse_gene():
    gene_cell = "GENE:G6PD"
    obs = parse_gene(gene_cell)
    exp = "G6PD"
    assert obs == exp


