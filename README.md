# Parser of PharmGKB's Allele Definition Tables
### example output 
* allele_definitions.tsv

|name|gene|chromosome|hgvs|start|end|rsid|variant_type|type|
|---|---|---|---|---|---|---|---|---|
|Reference|CACNA1S|chr1|['g.201091993G>A', 'g.201060815C>T']|['201091993', '201060815']|['201091993', '201060815']|['rs772226819', 'rs1800559']|['SNP', 'SNP']|['G', 'C']|
|c.520C>T|CACNA1S|chr1|['g.201091993G>A', 'g.201060815C>T']|['201091993', '201060815']|['201091993', '201060815']|['rs772226819', 'rs1800559']|['SNP', 'SNP']|['A', 'C']|
|c.3257G>A|CACNA1S|chr1|['g.201091993G>A', 'g.201060815C>T']|['201091993', '201060815']|['201091993', '201060815']|['rs772226819', 'rs1800559']|['SNP', 'SNP']|['G', 'T']
|...|...|...|...|...|...|...|...|...|
|rs9923231 variant (T)|VKORC1|chr16|['g.31096368C>T']|['31096368']|['31096368']|['rs9923231']|['SNP']|['T']|

* allele_definitions_name_hgvs_relation.tsv

|gene|name|hgvs|
|---|---|---|
|CACNA1S|c.520C>T|['g.201091993G>A']|
|CACNA1S|c.3257G>A|['g.201060815C>T']|
|...|...|...|
|CFTR|E831X|['g.117530975G>A', 'g.117594930G>T']|
|...|...|...|
|VKORC1|rs9923231 variant (T)|['g.31096368C>T']|

* allele_definitions_hgvs_name_relation.tsv

|gene|hgvs|name|
|---|---|---|
|CACNA1S|g.201091993G>A|['c.520C>T']|
|CACNA1S|g.201060815C>T|['c.3257G>A']|
|...|...|...|
|CFTR|g.117530975G>A|['E831X', 'R117H']|
|...|...|...|
|VKORC1|g.31096368C>T|['rs9923231 variant (T)']|

# create and activate conda environment
```shell
conda env create -f environment.yml
conda activate pgkb-dev
```

# install pgkb_hap_parser package
```shell
pip install -e .
```

# run pgkb_hap_parser package
```shell
pgkb_hap_parser parse -d DELIMITER_PATH -a ALLELE_DEFINITION_TABLE_PATH -o OUTPUT_PATH
```

# run test pgkb_hap_parser package
```shell
pytest-watch -c
```

# export conda environment
```shell
conda env export --no-builds | grep -v "prefix" > environment.yml
```
