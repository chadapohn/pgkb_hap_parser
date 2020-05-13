#!/bin/bash

# check if the conda environment exists
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pgkb-dev
if [ $? -eq 0 ]; then
   echo "conda environment 'pgkb-dev' does exist."
   echo "---> export conda environment."
   conda env export --no-builds | grep -v "prefix" > environment.yml
fi
