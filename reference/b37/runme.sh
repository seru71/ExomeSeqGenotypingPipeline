#!/bin/bash
#
# Recipe for reproducing interval-files for gene coverage calculation in GATK.
#
# Pawel Sztromwasser
# Feb, 2016
#


# Downloaded CCDS table from UCSC TableBrowser (hg b37), joined with name2 column from refGene (i.e. HGNC gene symbol).
# The downloaded table is in ccds.b37.jan2016

#
# 1. Create genes file (one record per CCDS, i.e. can be multiple records per gene symbol)
#
head -n1 ccds.b37.jan2016 | cut -f1-16 > ccds_genes.GRCh37.gatk
grep -v ^# ccds.b37.jan2016 | \
        sort -n -k5,5 | \
        sed 's/\tchr/\t/g' | \
        perl sortbyref.pl -k 3 - human_g1k_v37.fasta.fai | \
	sed 's/,$//g' | \
	awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$17"_"$2"\t"$14"\t"$15"\t"$16}' \
>> ccds_genes.GRCh37.gatk

#
# 2. Create exons file (one record per each exon-CCDS pair)
#
head -n1 ccds.b37.jan2016 | cut -f1-16 > ccds_exons.GRCh37.gatk
grep -v ^# ccds.b37.jan2016 | \
	sort -n -k5,5 | \
	sed 's/\tchr/\t/g' | \
	perl sortbyref.pl -k 3 - human_g1k_v37.fasta.fai | \
	sed 's/,$//g' | \
	python convert_ucsc_table_to_gatk_exons_table.py \
>> ccds_exons.GRCh37.gatk

#
# 3. Create BED file based on exons (no merging)
#
awk '{print $3"\t"$10"\t"$11"\t"$13"\t0\t"$4}' ccds_exons.GRCh37.gatk | \
	grep -v hg19.ccdsGene.chrom | \
	sed 's/,\t/\t/g' | \
	bedtools sort -i - \
> ccds_exons.not_merged.GRCh37.gatk.bed
#
# for merged intervals, pipe the above output through
#	bedtools merge -nms -i -

#
# 4. Chop the exons, so that in interval-based processing in GATK each overlapping fragment is attributed only to genes it belongs to. 
#    For details see chop_exons.py
#
python chop_exons.py ccds_exons.not_merged.GRCh37.gatk.bed > ccds_exons.chopped.GRCh37.gatk.bed

# 
# To confirm that intervals before and after chopping cover the same regions:
#
# $> bedtools intersect -a ccds_exons.chopped.GRCh37.gatk.bed -b ccds_exons.not_merged.GRCh37.gatk.bed -v
# $> bedtools intersect -b ccds_exons.chopped.GRCh37.gatk.bed -a ccds_exons.not_merged.GRCh37.gatk.bed -v
#
# No output means, that they don't differ.
#

