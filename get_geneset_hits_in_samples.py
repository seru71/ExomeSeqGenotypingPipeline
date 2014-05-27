#!/usr/bin/env python
"""
For a list of genes in GENE_SET_FILE and a list of files with variant calls 
(one sample per file; currently annovar avinput.hg19_multianno.with_omim.csv files) 
print out which genes from the gene set are found altered (per sample)

author: Pawel Sztromwasser (pawel.sztromwasser at k2.uib.no)  
"""

def get_genes(gene_set_file):
	""" Get the set of genes in interest """ 
	f = open(gene_set_file)
	genes = [e.strip() for e in f.xreadlines()]
	f.close()
	return set(genes)

from utility_functions import parenthesis_aware_split

def find_geneset_hits_in_samples(sample_files, geneset, gene_column_name = 'Gene.refGene'):
	"""
	Iterate over the per-sample variant tables (csv) and look for genes from the geneset.
	Returns a dict in form {sample: {gene1:[variant_line], gene2:[variant_line1, variant_line2]}}, where gene1 and gene2 belong to the geneset, 
	and variant_line is the entire line from csv file.
	As the same gene can appear twice with distinct variants, the dict[sample][gene] is a list.
	"""
	
	map={s:{} for s in sample_files}
	
	for fname in sample_files:
		sys.stderr.write('Processing '+ fname + '...')
		f = open(fname)
		
		header = quote_aware_split(table_in.readline().strip())
		gene_col_index = header.index(gene_column_name)
		
		for l in f.xreadlines():
			gene_entry = quote_aware_split(l)[gene_col_index]
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
					try:
						map[fname][gene] += [l]
					except KeyError:
						map[fname][gene] = [l]
										
		f.close()
        sys.stderr.write('done\n')
		
	return map



import sys
import glob

if __name__ == '__main__':
	
	if len(sys.argv) < 3:
		print 'Usage:', sys.argv[0], 'CSV_FILE(S) GENE_SET_FILE'
		print 'Example:', sys.argv[0],'sample_1234.avinput.hg19_multianno.with_omim.csv interesting_genes.txt'
		print 'Example:', sys.argv[0],'"annovar-files/*.avinput.hg19_multianno.with_omim.csv" interesting_genes.txt'
		sys.exit(0)
	
	variant_files = glob.glob(sys.argv[1])
	geneset = get_genes(sys.argv[2])
	
	map = find_geneset_hits_in_samples(variant_files, geneset, gene_column=1)

#	print '---------------'	
	for sample in sorted(map.keys()):
		print sample
		for gene in sorted(map[sample].keys()):
			print ''.join(map[sample][gene])
		print '------------------'
		
	
	
