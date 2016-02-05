#!/usr/bin/env python
#
# Usage: 
# Download refGene table from UCSC TableBrowser. 
# If you would like to use CCDS coordinates (one per gene), download CCDS->ccdsGene table joined with name2 column from refGene).
# The script is configured for CCDS+name2 table, and it will work with plain refGene only when modified.
# Note that some genes will have several transcripts/CCDs per gene.
#
# Example:
# $> head -n1 ccds.b37 > ccds_exons.sorted.b37
# $> grep -v ^# ccds.b37 | sort -n -k5,5 | perl sortbyref.pl -k 3 - human_genome.b37.fasta.fai | python convert_refseq_genes_to_exons.py >> ccds_exons.sorted.b37
#
# by Pawel Sztromwasser, 2016
#


import sys

# column indexes, 0-based
cols={'name':1, 'chr':2, 'exonCnt':8, 'starts':9, 'ends':10, 'score':11, 'name2':12, 'gene':16}

header='#bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        transcript      name2   cdsStartStat    cdsEndStat      exonFrames'

for l in sys.stdin.xreadlines():
        if l.find('#') == 0:
                print header
                continue
        lsplit = l.split('\t')

        # skip irregular chromosomes
        if lsplit[cols['chr']].find('_')>=0:
                continue
        starts = lsplit[cols['starts']].split(',')
	ends   = lsplit[cols['ends']].split(',')
	gene   = lsplit[cols['gene']].strip()
        for exon in range(0, int(lsplit[cols['exonCnt']])):
                print '\t'.join([lsplit[0],
                                gene,
                                '\t'.join(lsplit[cols['name']+1:cols['exonCnt']]),
                                '1',
                                starts[exon] + ',',
                                ends[exon] + ',',
                                lsplit[cols['name']],
                                gene+'_'+lsplit[cols['name']]+'_exon'+str(exon+1),
                                '\t'.join(lsplit[cols['name2']+1:-2]),
                                lsplit[-2].split(',')[exon]+','])

