from pgkb_hap_parser.parser import parse_deletion, parse_insertion, parse_snp

def test_parse_deletion():
	single = "g.94949283delA" # CYP2C9
	assert parse_deletion(single) == ("94949282", "94949283")

	multiple = "g.94942213_94942222delAGAAATGGAA" # CYP2C9
	assert parse_deletion(multiple) == ("94942212", "94942222")

def test_parse_insertion():
	single = "g.48041103_48041104insG" # NUDT15
	assert parse_insertion(single) == ("48041103", "48041104")

	multiple = "g.48037801_48037802insGAGTCG" # NUDT15
	assert parse_insertion(multiple) == ("48037801", "48037802")

def test_parse_snp():
	# Diallelic SNPs
	assert parse_snp("g.48037798G>A") == ("48037797", "48037798") #NUDT15

	# Multiallelic SNPs
	assert parse_snp("g.41016810C>A/T") == ("41016809", "41016810") # CYP2B6
	