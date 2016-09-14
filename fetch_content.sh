#!/bin/bash 

## 
# fetch markdown files and corresponding assets created by Rmarkdown 
# Args:
#   $1: path to directory

filename=$1
pagepath="pages/bioqc/"

read -p "Fetching \"$filename\" including assets. Do you want to continue? [y/n]" -n 1 -r
echo 
if [[ $REPLY =~ ^[Yy]$ ]]; then
    cp -vf $filename $pagepath
    assets="${filename/\.md/_files}"
    if [ -d $assets ]; then
        cp -Rfv $assets $pagepath
    else
        echo "No assets found." 
    fi
fi

##
# Adapt yaml header for jekyll
#
basename=$(basename $filename)
pagename=$pagepath$basename
head -n2 $pagename > $pagename.tmp
echo "permalink: ${basename/\.md/.html}" >> $pagename.tmp
tail -n+3 $pagename >> $pagename.tmp
mv $pagename.tmp $pagename
