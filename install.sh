#!/bin/bash

# check if the environment.yml exists
if [ -e "environment.yml" ]; then
    env=$(head -n 1 environment.yml | cut -f2 -d " ")
    # check if the conda environment exists
    source ~/anaconda3/etc/profile.d/conda.sh
    conda activate $env
    if [ $? -eq 0 ]; then
        echo "conda environment '$env' does exist."
        echo "---> activate and update conda environment."
        conda activate $env
        conda env update -f environment.yml -q
        echo "---> update pgkb_hap_parser package."
        pip install -e . -q
    else
        echo "conda environment '$env' doesn't exist."
        echo "---> create and activate conda environment."
        conda env create -f environment.yml -q
        conda activate $env
        echo "---> install pgkb_hap_parser package."
        pip install -e . -q
    fi
fi
