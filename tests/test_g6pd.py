import pandas as pd

from pgkb_hap_parser.parser import parse_gene

def test_parse_gene():
    gene_cell = "GENE:G6PD"
    obs = parse_gene(gene_cell)
    exp = "G6PD"
    assert obs == exp

def test_one_del_variant():
    # DPYD
    var_cells = {'chrom_hgvs_name': ["g.97450066delG"]}
    var_df = pd.DataFrame(data=var_cells)
    print(var_df)


# def test_insertion()

# def test_deletion()