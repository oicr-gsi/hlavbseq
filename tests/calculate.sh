#!/bin/bash
cd $1

# We have three output txt files, files with raw and formatted predictions

echo ".txt files:"
for t in *.txt;do md5sum $t;done 
