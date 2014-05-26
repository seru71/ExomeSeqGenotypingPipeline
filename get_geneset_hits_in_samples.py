#!/usr/bin/env python
"""
For a list of genes in GENE_SET_FILE and a list of files with variant calls 
(one sample per file; currently annovar avinput.refGene.variant_function files) 
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

def find_geneset_hits_in_samples(sample_files, geneset, gene_column = 1, delim='\t'):
	"""
	Iterate over the per-sample variant files and look for genes from the geneset.
	Returns a dict in form {sample: [gene1,gene2]}, where gene1 and gene2 belong to the geneset.
	The same gene appearing twice in the sample's list indicates two distinct variants in the gene.
	"""
	
	map={f:[] for f in sample_files}
	
	for fname in sample_files:
		sys.stderr.write('Processing '+ fname + '...')
		f = open(fname)
		for l in f.xreadlines():
			gene_entry = l.split('\t')[gene_column]
			
			gene_entry = gene_entry.strip().strip('"')		# clean the gene name
			genes=[gene_entry]
			if gene_entry.find(',') >= 0: 
				genes = parenthesis_aware_split(gene_entry,delim=',')	# split multi gene entries
			genes = [parenthesis_aware_split(gene,delim=';') for gene in genes]
			genes = set([gene for sublist in genes for gene in sublist])  # get unique gene ids only
			
			for gene in genes:
				if gene.find('(') > 0: 
					gene = gene[:gene.find('(')] # strip the transcript change in parenthesis
			
				if gene in geneset:
					map[fname] += [gene]
										
		f.close()
                sys.stderr.write('done\n')

		
	return map



import sys
import glob

if __name__ == '__main__':
	
	if len(sys.argv) < 3:
		print 'Usage:', sys.argv[0], 'AVINPUT_FILE(S) GENE_SET_FILE'
		print 'Example:', sys.argv[0],'sample_1234.avinput.refGene.variant_function interesting_genes.txt'
		print 'Example:', sys.argv[0],'"annovar-files/*.avinput.refGene.variant_function" interesting_genes.txt'
		sys.exit(0)
	
	variant_files = glob.glob(sys.argv[1])
	geneset = get_genes(sys.argv[2])
	
	map = find_geneset_hits_in_samples(variant_files, geneset, gene_column=1)

#	print '---------------'	
	for sample in sorted(map.keys()):
		print sample, ':', ','.join(map[sample])
	
	
