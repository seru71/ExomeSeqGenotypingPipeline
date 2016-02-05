
Gene-coverage stats in GATK
---------------------------


The set of files used for gene coverage calculation can be recreated by running:
./runme.sh

The ccds_genes.GRCh37.gatk and ccds_exons.GRCh.gatk files are meant to be provided in -geneList arg for GATK's DepthOfCoverage.
To get the gene- or/and exon-based depth computed properly, use " -L ccds_exons.chopped.GRCh37.gatk.bed --interval_merging OVERLAPPING_ONLY ".
Example:
	java -jar GenomeAnalysisTK.jar -R $REFERNCE -T DepthOfCoverage -I $BAM /
		-geneList ccds_genes.GRCh37.gatk -L ccds_exons.chopped.GRCh37.gatk.bed --interval_merging OVERLAPPING_ONLY /
		--omitDepthOutputAtEachBase --omitLocusTable -ct 5 -ct 10 -ct 20 /
		-o $BAM

Also, to get proper stats for overlapping genes, make sure that you are using GATK with this patch:
https://github.com/broadgsa/gatk/pull/16
https://github.com/seru71/gatk/commit/d8f4c4667b369d1f8f8502a157adad0f87bee563


--------------------

Paweł Sztromwasser
Feb 2016
