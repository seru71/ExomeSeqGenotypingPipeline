#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher, modified by Paweł Sztromwasser
:Contact: mkircher [at] uw.edu, pawel.sztromwasser [at] ii.uib.no
:Date: *23.03.2014, modified 30.04.2014

:Description:

This script reads VCF from standard input and places the retrieved CADD 
scores in the vcf INFO field, spitting the vcf to standard output.
---------- this is not happening right now --------------
while generating a file with variants that could not be found in the downloaded file. Use it with a 
gzip compressed VCF and a downloaded score file (incl. downloaded index) as follows:

gunzip -c input.vcf.gz | python extractScoresVCF.py -p CADDscores.tsv.gz | gzip -c > scores.tsv.gz

Or if you gave execute permission to the script (chmod +x extractScoresVCF.py):

gunzip -c input.vcf.gz | ./extractScoresVCF.py -p CADDscores.tsv.gz | gzip -c > scores.tsv.gz

In case your input.vcf is not compressed by gzip replace 
"gunzip -c input.vcf.gz" by "cat input.vcf".

It will create an uncompressed indels_out.vcf in the same folder
which can then be submitted to the CADD webserver for scoring 
the remaining variants.

"""

import sys, os
import pysam
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p","--path", dest="path", help="Path to scored genome (default './whole_genome_SNVs_anno.tsv.gz')",default="./whole_genome_SNVs_anno.tsv.gz")
parser.add_option("--indels_out", dest="indels_out", help="Write indels and other not found variants to file (default 'indels_out.vcf')",default="indels_out.vcf")
(options, args) = parser.parse_args()

# VCF FIELDS
fchr = 0
fpos = 1
fdbsnp = 2
fref_allele = 3
falt_allele = 4
fgeno_qual = 5
fflag = 6
finfo = 7
fformat = 8
fvalues = 9

# included info line in the header
cadd_info_line = '##INFO=<ID=CADD,Number=A,Type=Float,Description="CADD score of variant deleteriousness; two values per alt allele separated by a comma">'

if os.path.exists(options.path) and os.path.exists(options.path+".tbi"):
  filename = options.path
  sys.stderr.write("Opening %s...\n"%(filename))
  regionTabix = pysam.Tabixfile(filename,'r')
else:
  sys.stderr.write("Require valid file with scored variants.\n")
  sys.exit()

found_first_info=False

for line in sys.stdin:
  if line.startswith('#'): 
	if line.startswith('##INFO='): found_first_info=True
	elif found_first_info:
		print cadd_info_line
		found_first_info=False
	print line,
	continue

  fields = line.rstrip().split('\t')
  chrom = fields[fchr]
  pos = int(fields[fpos])
  lref = fields[fref_allele]
  present_alleles = set(fields[falt_allele].split(','))

  sys.stdout.write('\t'.join(fields[:finfo+1]) + ';CADD=')

  first = True 
  for allele in present_alleles:
    found = False
    for VEPline in regionTabix.fetch(chrom,pos-1,pos):
      vfields = VEPline.rstrip().split('\t')
      if len(vfields) < 6: continue
      
      if vfields[2] == lref and vfields[3] == allele:
	if not first: sys.stdout.write(',')
        sys.stdout.write(vfields[5]+'('+vfields[4]+')')
        found = True
	break

    if not found:
    	if not first: sys.stdout.write(',')
      	sys.stdout.write('NA')

    first = False

  sys.stdout.write('\t'+'\t'.join(fields[finfo+1:])+'\n')
