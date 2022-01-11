#!/bin/bash
cd $1

# We have three output txt files, files with raw and formatted predictions

echo ".txt file:"
for c in *txt;do md5sum $c;done 
