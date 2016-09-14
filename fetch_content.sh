#!/bin/bash 

## 
# fetch markdown files and corresponding assets created by Rmarkdown 
# Args:
#   $1: path to directory

filename=$1

read -p "Fetching \"$filename\" including assets. Do you want to continue? [y/n]" -n 1 -r
echo 
if [[ $REPLY =~ ^[Yy]$ ]]; then
    cp -vf $filename ./pages/bioqc
    assets="${filename/\.md/_files}"
    if [ -d $assets ]; then
        cp -Rfv assets ./pages/bioqc
    else
        echo "No assets found." 
    fi
fi
