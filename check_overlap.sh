#!/bin/bash
#
# Check overlap between any given number of tab-separated variant tables. 
# Two records are considered unique if they have identical first 8 columns (tunnable by UNIQUR_COLS variable)
# 
# Usage: check_overlap.sh sample1.tsv sample2.tsv sample3.tsv
#
# Result is a table of variants that is shared by at least 2 samples (tables). 
# The first column (#Samples) shows the number of samples the variant was found in. The rest is the first 8 columns of the tables (indetical for #Samples of samples)
#
#
UNIQUE_COLS=8
cat $* | cut -f1-${UNIQUE_COLS} | sort -n | uniq -cd | sed 's/^ \+//g' | sed 's/[[:digit:]]\+ Chr\t/#Samples Chr\t/g' | sed 's/ /\t/g'
