
# ExomeSeq Genotyping Pipeline 



## Pipeline

The pipeline consists of SAMtools, Picard-tools, GATK, and VCFTools 



## Setup

### Pipeline

1. Install Python if you don't have it from before, and a cool python package - `Ruffus` (http://www.ruffus.org.uk/). 
You might also need following packages: optparse, logging, shutil

2. Clone the pipeline repository:
`git clone https://github.com/seru71/ExomeSeqGenotypingPipeline.git <PIPELINE_HOME>`

The pipeline is ready now, but you will need all of its components to perform the analysis.

### Pipeline components

Install following tools:
1. SAMtools (https://github.com/samtools/samtools)
2. Picard-tools (http://broadinstitute.github.io/picard)
3. GATK (https://www.broadinstitute.org/gatk/download)
4. VCFtools (https://vcftools.github.io)

Optional tool for some QC tasks:
5. QualiMap (http://qualimap.bioinfo.cipf.es)

### Reference data

1. Download resource bundle from GATK (1000Genomes version preferably):
https://www.broadinstitute.org/gatk/download

2. Unpack it in a directory of your choice (e.g. reference), and prepare the reference genome for use with GATK as recommended:
http://gatkforums.broadinstitute.org/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference
 
3. Download/create BED files representing target regions of you exome capture assay. The pipeline supports two interval files:
`capture` - used for QC metrics (QualiMap, depth of coverage reports, etc.)
`exome`   - usually superset of `capture`; used to limit data processing to these regions. No variants outside of these intervals will be reported.
They can be the exact same file, but sometimes one wants to add flanking sequence to the capture (e.g. +/-10bp) when filtering variants. 
Both files should be prepared in advance.     



## Usage

*   Running the script

    You can run the script using `python ~/tools/ExomeSeqGenotypingPipeline/pipeline_multisample.py`.
    A list of possible options will be presented. The only required option is `--pipeline_settings`,
    which accepts a config file with paths to input bam files, resources (e.g. reference genome, database, capture regions), 
    and right versions of software. See an exemplary file for all required options 
    in ~/tools/ExomeSeqGenotypingPipeline/pipeline_settings.cfg 
  
    If you want to follow progress of the script, use the verbose option (`-vvvvv`).
    In order to use multithreading, use the `-j` option (`-j 12`).

*   Outputs

    The script will create one or more directories in the current folder, named after
    the bam file prefix (`xxx.bam` will create directory `xxx`).

    After finishing a few files will be in that directory, including the cleaned up 
    bam file (`xxx.gatk.bam`), a .gvcf file, some quality control files (e.g. coverage statistics),
    and the vcf file with variants restricted to the exome.

    In the directory where the script is run, additional files will be created with the
    multisample results (`multisample.gatk.vcf` for the multisample snp and indels,
    and `multisample.gatk.analysisReady.exome.vcf` for the variants restricted to the exome bed-file).

*   Typical usage

    For running the genotyping analysis using 12 cpus (half of astrakan).

	mkdir ~/results/new_project
	cd ~/results/new_project
	cp ~/tools/ExomeSeqGenotypigPipeline/pipeline_settings.cfg .
	vi pipeline_settings.cfg  ## insert correct paths to the bam files, ref genome, executables, etc.
	pipeline_multisample.py -s pipeline_settings.cfg -t split_snps -vvvvv -j 12


### pipeline_annovar.py

The script for annotating variants with location (gene, exon, intron), effect on protein (synonymous, nonsynonymous, stopgain), 
and various pathogenicity scores (OMIM, PolyPhen, CADD). 






