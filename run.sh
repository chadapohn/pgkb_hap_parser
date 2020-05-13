#!/bin/bash

# check if the conda environment exists
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pgkb-dev
if [ $? -eq 0 ]; then
    echo "conda environment 'pgkb-dev' does exist."
    echo "---> run pgkb_hap_parser package."
    delimiter="$PWD/data/allele_definition_table/delimiter.tsv"
    allele_definition_table="$PWD/data/allele_definition_table"
    output="$PWD/out"
    pgkb_hap_parser parse -d $delimiter -a $allele_definition_table -o $output
fi
