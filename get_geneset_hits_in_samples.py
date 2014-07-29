#!/usr/bin/env python
"""
For a list of genes in GENE_SET_FILE and a list of files with variant calls 
(one sample per file; currently annovar avinput.hg19_multianno.with_omim.csv files and .variant_function files are supported)
print out which genes from the gene set are found altered (per sample)

author: Pawel Sztromwasser (pawel.sztromwasser at k2.uib.no)  
"""

def get_genes(gene_set_file):
	""" Get the set of genes in interest """ 
	f = open(gene_set_file)
	genes = [e.strip() for e in f.xreadlines()]
	f.close()
	return set(genes)

from utility_functions import quote_aware_split, parenthesis_aware_split

def find_geneset_hits_in_samples(sample_files, geneset, gene_column_name = 'Gene.refGene', delim=',', has_header=True, fields=None):
	"""
	Iterate over the per-sample variant tables (csv, tsv) and look for genes from the geneset.
	The gene column is given in the gene_column_name argument: either by name (if has_header==True), or by index of the column (if has_header==False).
	Returns a dict in form {sample: {gene1:[variant_line], gene2:[variant_line1, variant_line2]}}, where gene1 and gene2 belong to the geneset, 
	and variant_line is the entire line from csv file.
	As the same gene can appear twice with distinct variants, the dict[sample][gene] is a list.
	"""
	
	map={s:{} for s in sample_files}
	
	for fname in sample_files:
		sys.stderr.write('Processing '+ fname + '...')
		f = open(fname)
		
		gene_col_index = -1
		if has_header:
			header = quote_aware_split(f.readline().strip(), delim)
			gene_col_index = header.index(gene_column_name)
		else:
			gene_col_index = int(gene_column_name)
			
		
		for l in f.xreadlines():
			gene_entry = quote_aware_split(l, delim)[gene_col_index]
			gene_entry = gene_entry.strip().strip('"') 	# clean the gene name		
			genes=[gene_entry]						
			if gene_entry.find(',') >= 0: 
				genes = parenthesis_aware_split(gene_entry,delim=',')	# split multi gene entries
			genes = [parenthesis_aware_split(gene,delim=';') for gene in genes]
			genes = set([gene for sublist in genes for gene in sublist])  # get unique gene ids only
			
			for gene in genes:
				if gene.find('(') > 0: 
					gene = gene[:gene.find('(')] # strip the transcript change in parenthesis
			
				if gene in geneset:
					# clean the records
					variant_record = '\t'.join([e.strip().strip('"') for e in quote_aware_split(l, delim)])
					# if requested, select a subset of fields
					if fields != None:
						variant_record = '\t'.join([variant_record.split('\t')[i] for i in fields])

					try:
						map[fname][gene] += [variant_record]
					except KeyError:
						map[fname][gene] = [variant_record]
										
		f.close()
	        sys.stderr.write('done\n')
		
	return map



import sys
import glob
from os import path

if __name__ == '__main__':
	
	if len(sys.argv) < 3:
		print 'Usage:', sys.argv[0], 'CSV_FILE(S) GENE_SET_FILE'
		print 'Example:', sys.argv[0],'sample_1234.avinput.hg19_multianno.with_omim.csv interesting_genes.txt'
		print 'Example:', sys.argv[0],'"annovar-files/*.avinput.hg19_multianno.with_omim.csv" interesting_genes.txt'
                print 'Example:', sys.argv[0],'"annovar-files/*.avinput.variant_function" interesting_genes.txt'
		sys.exit(0)
	
	variant_files = glob.glob(sys.argv[1])
	geneset = get_genes(sys.argv[2])
	
	# settings for cnv files
	delim = ','
	has_header = True
	gene_column_name = 'Gene.refGene'
	selected_columns=[19] + range(0,19) + [27,29,31]
        header=['Sample','Gene','Zygozity','Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','ExonicFunc.refGene','AAChange.refGene',
                '1000g2012apr_eur','1000g2012apr_amr','1000g2012apr_asn','1000g2012apr_afr','dbSNP138','OMIM_phenotype','SIFT_avsift',
                'clinvar_20140211','PolyPhen2_ljb23_pp2hvar','CADDgt10','QUAL','vcf_INFO','GT:AD:DP:GQ:PL']

	# settings for annovar variant_function files
	if variant_files[0].split('.')[-1] != 'csv':
		delim = '\t'
		has_header = False
		gene_column_name = '1'
		selected_columns = range(0,10)+[12,17,19]
	        header=['Sample','Gene','Func.refGene','Gene.refGene','Chr','Start','End','Ref','Alt','Zygozity','QUAL','depth','dbSNP','vcf_INFO','GT:AD:DP:GQ:PL']


	map = find_geneset_hits_in_samples(variant_files, geneset, gene_column_name, delim, has_header, selected_columns)
	"""
	print '\t'.join(header)
	for sample in sorted(map.keys()):
		for gene in sorted(map[sample].keys()):
			sample_gene = path.basename(sample).split('.')[0]+'\t'+gene+'\t'
			print sample_gene + ('\n'+sample_gene).join(map[sample][gene])
		
	"""

	for sample in sorted(map.keys()):
#		variant_count=0
#		genes_str=''
#                for gene in map[sample].keys():
#			variant_count+=len(map[sample][gene])
#                        genes_str+=';'+gene+'('+str(len(map[sample][gene]))+')'
		variant_count = sum([len(map[sample][gene]) for gene in map[sample].keys()])
		genes_str = ';'.join([gene+'('+str(len(map[sample][gene]))+')' for gene in map[sample].keys()])
		
		print '\t'.join([path.basename(sample), str(variant_count), genes_str])		
