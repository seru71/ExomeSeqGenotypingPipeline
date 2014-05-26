#!/usr/bin/env python
"""
For a list of genes in GENE_SET_FILE and a list of files with variant calls 
(one sample per file; currently annovar avinput.refGene.variant_function files) 
print out which genes from the gene set are found altered (per sample)

author: Paweł Sztromwasser (pawel.sztromwasser@k2.uib.no  
"""

def getGenes(gene_set_file):
	""" Get the set of genes in interest """ 
	f = open(gene_set_file)
	genes = [e.strip() for e in f.xreadlines()]
	f.close()
	return set(genes)


def parenthesis_aware_split(string, delim=',', open_par='(', close_par=')'):
	"""	Split outside of parenthesis (i.e. ignore delimiters within parenthesis."""
	out = []
	s = ''
	open_parenthesis=0
	for c in string:
		if c == open_par: 
			open_parenthesis+=1
		if c == close_par and open_parenthesis > 0:
			open_parenthesis-=1
		if c == delim and open_parenthesis==0:
			out += [s]
			s = ''
		else: 
			s += c
	return out + [s]


def findGenesInSamples(sample_files, geneset, gene_column = 1, delim='\t'):
	"""
	Iterate over the per-sample variant files and look for genes from the geneset.
	Returns a dict in form {sample: [gene1,gene2]}, where gene1 and gene2 belong to the geneset.
	The same gene appearing twice in the sample's list indicates two distinct variants in the gene.
	"""
	
	map={f:[] for f in sample_files}
	
	for fname in sample_files:
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
		
	return map



import sys
import glob

if __name__ == '__main__':
	
	if len(sys.argv) < 2:
		print 'Usage: getGenesetHitsInSample.py AVINPUT_FILE(S) GENE_SET_FILE'
		print 'Example: getGenesetHitsInSample.py sample_1234.avinput.refGene.variant_function interesting_genes.txt'
		print 'Example: getGenesetHitsInSample.py annovar-files/*.avinput.refGene.variant_function interesting_genes.txt'
	
	avinput_files = glob.glob(sys.argv[0])
	geneset = getGenes(sys.argv[1])
	
	map = mapGenesInSamples(avinput_files, geneset, gene_column=1)
	
	for sample in map:
		print sample, ':', ','.join(map[sample])
	
	