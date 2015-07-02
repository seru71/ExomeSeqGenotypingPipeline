#!/usr/bin/env python
"""

    pipeline_annovar.py
            [--bamdir PATH]
                        [--log_file PATH]
                        [--verbose]
                        [--target_tasks]
                        [--jobs]
                        [--just_print]
                        [--flowchart]
                        [--key_legend_in_graph]
                        [--forced_tasks]

"""
import sys
import os
import glob
import string


if __name__ == '__main__':
    from optparse import OptionParser
    import StringIO

    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --settings PIPELINE_SETTINGS.CFG --groups NUMBER [more_options]")
    parser.add_option("-s", "--settings", dest="pipeline_settings",
                        metavar="FILE",
                        type="string",
                        help="File containing all the settings for the analysis.")
                        


    #
    #   general options: verbosity / logging
    #
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="count", default=0,
                      help="Print more verbose messages for each additional verbose level.")
    parser.add_option("-L", "--log_file", dest="log_file",
                      metavar="FILE",
                      type="string",
                      help="Name and path of log file")




    #
    #   pipeline
    #
    parser.add_option("-t", "--target_tasks", dest="target_tasks",
                        action="append",
                        default = list(),
                        metavar="JOBNAME",
                        type="string",
                        help="Target task(s) of pipeline.")
    parser.add_option("-j", "--jobs", dest="jobs",
                        default=1,
                        metavar="N",
                        type="int",
                        help="Allow N jobs (commands) to run simultaneously.")
    parser.add_option("-n", "--just_print", dest="just_print",
                        action="store_true", default=False,
                        help="Don't actually run any commands; just print the pipeline.")
    parser.add_option("--flowchart", dest="flowchart",
                        metavar="FILE",
                        type="string",
                        help="Don't actually run any commands; just print the pipeline "
                             "as a flowchart.")

    #
    #   Less common pipeline options
    #
    parser.add_option("--key_legend_in_graph", dest="key_legend_in_graph",
                        action="store_true", default=False,
                        help="Print out legend and key for dependency graph.")
    parser.add_option("--forced_tasks", dest="forced_tasks",
                        action="append",
                        default = list(),
                        metavar="JOBNAME",
                        type="string",
                        help="Pipeline task(s) which will be included even if they are up to date.")
    parser.add_option("--rebuild_mode", dest="rebuild_mode",
                        action="store_false", default=True,
                        help="gnu_make_maximal_rebuild_mode")

    # get help string
    f =StringIO.StringIO()
    parser.print_help(f)
    helpstr = f.getvalue()
    (options, remaining_args) = parser.parse_args()


    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Change this if necessary                  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    #
    #   Add names of mandatory options,
    #       strings corresponding to the "dest" parameter
    #       in the options defined above
    #
    mandatory_options = ['pipeline_settings']

    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Change this if necessary                  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    def check_mandatory_options (options, mandatory_options, helpstr):
        """
        Check if specified mandatory options have been defined
        """
        missing_options = []
        for o in mandatory_options:
            if not getattr(options, o):
                missing_options.append("--" + o)
    
        if not len(missing_options):
            return
    
        raise Exception("Missing mandatory parameter%s: %s.\n\n%s\n\n" %
                        ("s" if len(missing_options) > 1 else "",
                         ", ".join(missing_options),
                         helpstr))
    check_mandatory_options (options, mandatory_options, helpstr)





    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Get pipeline settings from a config file  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    import ConfigParser

    config = ConfigParser.ConfigParser()
    config.read(options.pipeline_settings)
    # inputs 
    try: input_vcfs = config.get('Inputs','input-vcfs')
    except (ConfigParser.NoOptionError):
        sys.stderr.write('No input-vcfs setting in config file. Recreating vcf files path from bam files: ')
        prefix = os.path.splitext(os.path.basename(config.get('Inputs','input-bams')))[0]
        input_vcfs = os.path.join(prefix, prefix+'.exome.vcf') 
        sys.stderr.write(input_vcfs+'\n')


    # reference dbs
    annovar_human_db = config.get('Resources','annovar-humandb-dir')
    annovar_1000genomes_eur = config.get('Resources','annovar-1000genomes-eur')
    annovar_1000genomes_eur_maf = config.get('Resources','annovar-1000genomes-eur-MAF-cutoff')
    annovar_inhouse_db = config.get('Resources','annovar-inhouse-db')
    omim_gene_phenotype_map_file = config.get('Resources','omim_gene_phenotype_map')



#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   imports


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

from ruffus import *
import subprocess
import resource


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Logger


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

if __name__ == '__main__':
    import logging
    import logging.handlers

    MESSAGE = 15
    logging.addLevelName(MESSAGE, "MESSAGE")

    def setup_std_logging (logger, log_file, verbose):
        """
        set up logging using programme options
        """
        class debug_filter(logging.Filter):
            """
            Ignore INFO messages
            """
            def filter(self, record):
                return logging.INFO != record.levelno

        class NullHandler(logging.Handler):
            """
            for when there is no logging
            """
            def emit(self, record):
                pass

        # We are interesting in all messages
        logger.setLevel(logging.DEBUG)
        has_handler = False

        # log to file if that is specified
        if log_file:
            handler = logging.FileHandler(log_file, delay=False)
            handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)6s - %(message)s"))
            handler.setLevel(MESSAGE)
            logger.addHandler(handler)
            has_handler = True

        # log to stderr if verbose
        if verbose:
            stderrhandler = logging.StreamHandler(sys.stderr)
            stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
            stderrhandler.setLevel(logging.DEBUG)
            if log_file:
                stderrhandler.addFilter(debug_filter())
            logger.addHandler(stderrhandler)
            has_handler = True

        # no logging
        if not has_handler:
            logger.addHandler(NullHandler())


    #
    #   set up log
    #
    module_name = "exome"
    logger = logging.getLogger(module_name)
    setup_std_logging(logger, options.log_file, options.verbose)

    #
    #   Allow logging across Ruffus pipeline
    #
    def get_logger (logger_name, args):
        return logger

    from ruffus.proxy_logger import *
    (logger_proxy,
     logging_mutex) = make_shared_logger_and_proxy (get_logger,
                                                    module_name,
                                                    {})


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Functions


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


def setlimits():
    """Set maximum meomory to be used by a child process"""
    resource.setrlimit(resource.RLIMIT_AS, (100000000000,100000000000))

def run_cmd(cmd_str):
    """
    Throw exception if run command fails
    """
    process = subprocess.Popen(cmd_str, stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE, shell = True, preexec_fn=setlimits)
    stdout_str, stderr_str = process.communicate()
    if process.returncode != 0:
        raise Exception("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                            (cmd_str, stdout_str, stderr_str, process.returncode))


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Pipeline


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#       Put pipeline code here


def generate_parameters():
    vcfs = glob.glob(input_vcfs)
    for f in vcfs:
        # get the sample id as prefix
        prefix = os.path.basename(f)[0:-len('.exome.vcf')]
        yield [f, os.path.join('annotated-with-annovar', prefix+'.avinput')]


#@transform(vcf_files, formatter('.*/(?P<SAMPLE_ID>[^/]+).exome.vcf'),
#                        '{subpath[0][1]}/annotated-with-annovar/{SAMPLE_ID[0]}.avinput')
@files(generate_parameters)
def prepare_annovar_inputs(vcf, output):
    try: os.mkdir('annotated-with-annovar')
    except (OSError): pass # dir exists
 
    run_cmd("convert2annovar.pl {vcf} -format vcf4 -withzyg -includeinfo -outfile {out} \
        ".format(vcf=vcf, out=output))
    
   
@transform(prepare_annovar_inputs, suffix('.avinput'), 
                                        ['.avinput.hg19_EUR.sites.2012_04_filtered', 
                                         '.avinput.hg19_EUR.sites.2012_04_dropped'])
def filter_common_1000genomes(input, outputs):
    """ filter common 1000 genomes variants """
    run_cmd("annotate_variation.pl -build hg19 -filter -dbtype {eur1kg} \
        -maf {maf} -outfile {output_prefix} {input_file} {annodb}".format(
        eur1kg=annovar_1000genomes_eur,
        maf=annovar_1000genomes_eur_maf, 
        output_prefix=input,
        input_file=input, 
        annodb=annovar_human_db))

@transform(filter_common_1000genomes, suffix('.hg19_EUR.sites.2012_04_filtered'),
                                        ['.hg19_EUR.sites.2012_04_filtered.common_inhouse_filtered',
                                         '.hg19_EUR.sites.2012_04_filtered.common_inhouse_dropped'])
def filter_common_inhouse(inputs, outputs):
    filtered = inputs[0]                      # take only the filtered file, leave dropped

    """ filter variants found in the inhouse database. OBS! output specifies the filename after rename """    
    run_cmd("annotate_variation.pl -build hg19 -filter -dbtype generic -genericdbfile {inhouse} \
        -outfile {outfile} {input} {annodb}".format(
        inhouse=annovar_inhouse_db, 
        outfile=filtered, 
        input=filtered, 
        annodb=annovar_human_db))

    for output_file in outputs: 
        os.rename(output_file.replace('common_inhouse','hg19_generic'), output_file)



def get_stats_on_prefiltered_variants(input, outputs, cleanup=True):    

    run_cmd("annotate_variation.pl -buildver hg19 -outfile {outfile_prefix} {input_file} {annodb}".format(
        outfile_prefix=input, 
        input_file=input, 
        annodb=annovar_human_db))

    # calculate stats on files created by annovar - output files without ".stats" suffix
    run_cmd("cut -f 1 {f} | sort | uniq -c > {f}.stats".format(f=outputs[0][:-len('.stats')]))
    run_cmd("cut -f 2 {f} | sort | uniq -c > {f}.stats".format(f=outputs[1][:-len('.stats')]))
    # remove the annovar files
    if cleanup:
        os.remove(outputs[0][:-len('.stats')])
        os.remove(outputs[1][:-len('.stats')])


@transform(filter_common_inhouse, suffix('.common_inhouse_filtered'), 
                                           ['.common_inhouse_filtered.variant_function',
                                         '.common_inhouse_filtered.exonic_variant_function',
                                         '.common_inhouse_filtered.variant_function.stats',
                                         '.common_inhouse_filtered.exonic_variant_function.stats'])
def annotate_function_of_rare_variants(inputs, outputs):
    """ annotate functional change in rare variants """
    filtered = inputs[0]              # use only the filtered input file, leave dropped
    get_stats_on_prefiltered_variants(input=filtered, outputs=outputs[3:4], cleanup=False)


@transform(annotate_function_of_rare_variants, 
           formatter(".*/(?P<SAMPLE_ID>[^/]+).avinput.hg19_EUR.sites.2012_04_filtered.common_inhouse_filtered.variant_function", None, None, None),
           ['{path[0]}/annotated-tables/{SAMPLE_ID[0]}.rare_coding_and_splicing.avinput', 
            '{path[0]}/annotated-tables/{SAMPLE_ID[0]}.rare_coding_and_splicing.avinput.hg19_multianno.csv'])
def produce_variant_annotation_table(inputs, outputs):
    """ produce a table of various annotations per variant """
        
    dir = 'annotated-with-annovar/annotated-tables'
    try: os.mkdir(dir)
    except (OSError): pass # dir exists
    
    avinput = outputs[0]
    rare_var_fun = inputs[0]
    rare_ex_var_fun = inputs[1]
    
    # create input file for the table_annotation script
    f_out = open(avinput,'w')
    f = open(rare_ex_var_fun)
    for l in f.xreadlines():
        lsplit=l.split('\t')
        if lsplit[1] != 'synonymous SNV':
            f_out.write(string.join(lsplit[3:],'\t'))
    f.close()
    f = open(rare_var_fun)
    for l in f.xreadlines():
        lsplit=l.split('\t')
        if 'splicing' in (lsplit[0]).split(';'):
            f_out.write(string.join(lsplit[2:],'\t'))
    f.close()
    f_out.close()
    
    # annotate all variants selected above
    run_cmd("table_annovar.pl -protocol refGene,1000g2012apr_eur,1000g2012apr_amr,1000g2012apr_asn,1000g2012apr_afr,snp138,avsift,clinvar_20140211,ljb23_pp2hvar,caddgt10 \
            -operation g,f,f,f,f,f,f,f,f,f -arg \'-splicing 4\',,,,,,,,,\'-otherinfo\' -nastring NA -build hg19 -csvout -otherinfo \
            -outfile {output_prefix} {input} {db}".format(
        output_prefix=avinput, 
        input=avinput, 
        db=annovar_human_db))


#
# include omim phenotypes

def get_omim_gene_phenotype_map(omim_file):
    gene_col=6
    pht_col=12
    map_pht={}
    f = open(omim_file)
    for l in f.xreadlines():
        lsplit=l.split('|')
        
        # ignore lines with no phenotype
        if lsplit[pht_col-1].strip()=='':
            continue

        genes, phenotype = lsplit[gene_col-1], lsplit[pht_col-1]
        for gene in genes.split(','):
            gene = gene.strip()
            if gene == '': continue
            try:
                map_pht[gene] = map_pht[gene]+'|'+phenotype.strip()
            except KeyError:
                map_pht[gene] = phenotype.strip()
        
    f.close()
    return map_pht


omim_gene_phenotype_map = get_omim_gene_phenotype_map(omim_gene_phenotype_map_file)

from utility_functions import quote_aware_split, parenthesis_aware_split

@transform(produce_variant_annotation_table, formatter(), '{path[1]}/{basename[1]}.with_omim.csv')
def include_omim_phenotype_annotation(inputs, output_table, gene_column=7, omim_column=15, delim=','):
    """ include OMIM phenotype into the annotation table """
    table_in = open(inputs[1],'r')
    table_out = open(output_table,'w')

    # header
    header_in=quote_aware_split(table_in.readline(), delim)
    if omim_column <= 0:        
        omim_column=len(header_in)+1
    header_out = header_in[:omim_column-1] + ['omim_phenotype'] + header_in[omim_column-1:]
    table_out.write(delim.join(header_out))

    # the rest of the table
    for l in table_in.xreadlines():
        lsplit = quote_aware_split(l,delim)
        gene = lsplit[gene_column-1].strip('"')

        # the gene record can be a list (e.g. overlapping genes), so it needs to be split
        genes = [gene]
        if gene.find(',')>=0:
            genes = parenthesis_aware_split(gene, delim=',')
        genes = [parenthesis_aware_split(gene,delim=';') for gene in genes]
        genes = set([gene for sublist in genes for gene in sublist])  # unlist and get unique gene ids only
        
        for gene in genes:
            # if present, remove suffix in parenthesis
            if gene.find('(') >= 0: 
                gene = gene[:gene.find('(')]
            
            # put the variant record in the map
            try:
                omim_phenotype = omim_gene_phenotype_map[gene]
            except KeyError:
                omim_phenotype = 'NA'

        table_out.write(delim.join(lsplit[:omim_column-1] + ['"'+omim_phenotype+'"'] + lsplit[omim_column-1:]) )

    table_in.close()
    table_out.close()

@transform(include_omim_phenotype_annotation, formatter(), '{path[0]}/{basename[0]}.recessive.csv')
def extract_recessive_disorder_candidates(input, output, gene_column_name='Gene.refGene', zygozity_column_name='Otherinfo', delim=','):
    """ extract a part of the annotated table that contains candidates for recessive disorders """ 
    table_in = open(input,'r')
    
    # get right column indexes
    header = quote_aware_split(table_in.readline().strip(), delim)
    gene_col_index     = header.index(gene_column_name)
    zygozity_col_index = header.index(zygozity_column_name) 
    
    variant_records_per_gene={}
    for l in table_in.xreadlines():
        lsplit = quote_aware_split(l,delim)
        
        gene = lsplit[gene_col_index].strip('"')
        
        # the gene record can be a list (e.g. overlapping genes), so it needs to be split
        genes = [gene]
        if gene.find(',')>=0:
            genes = parenthesis_aware_split(gene, delim=',')
        genes = [parenthesis_aware_split(gene,delim=';') for gene in genes]
        genes = set([gene for sublist in genes for gene in sublist])  # unlist and get unique gene ids only
        
        for gene in genes:
            # if present, remove suffix in parenthesis
            if gene.find('(') >= 0: gene = gene[:gene.find('(')]
            # put the variant record in the map
            try:
                variant_records_per_gene[gene] += [l]
            except KeyError: 
                variant_records_per_gene[gene] = [l]
    
    table_in.close()
    
    # write the table with candidates for recessive inheritance model
    table_out = open(output,'w')
    table_out.write(delim.join(header)+'\n')
    
    # iterate over the genes and select...
    for gene in variant_records_per_gene:
        # ...these with 2 or more variants...
        if len(variant_records_per_gene[gene]) >= 2:
            for l in variant_records_per_gene[gene]: 
                table_out.write(l)
        # or homozygous variants
        else:
            lsplit = quote_aware_split(variant_records_per_gene[gene][0])
            if lsplit[zygozity_col_index].find('"hom\t') >= 0:
                table_out.write(variant_records_per_gene[gene][0])
    
    table_out.close()
    
    

#
# QC on variant level
##

@merge(prepare_annovar_inputs, ['hetz_per_chr.tsv','homz_per_chr.tsv'])
def count_hetz_and_homz_per_chr(infiles, table_files):
    """ produce a table of sample vs chromosome counts of hetero- and homozygotic variants """

    # the variant calls are expected to appear sorted by chromosome name in the following order
    # any other order should trigger an exception
    chromosomes = ['1','2','3','4','5','6','7','8','9',
                   '10','11','12','13','14','15','16','17','18','19',
                   '20','21','22','X','Y','MT',
                   'GL000191.1', 'GL000228.1', 'GL000209.1',  
                   'GL000223.1','GL000222.1', 'GL000194.1']

    hetz = open(table_files[0], 'w')
    homz = open(table_files[1], 'w')
    hetz.write('sample\t'+string.join(chromosomes,'\t')+'\n')
    homz.write('sample\t'+string.join(chromosomes,'\t')+'\n')
    for fname in infiles:
        sample_id=os.path.basename(fname).split('.')[0]
        hetz.write(sample_id)
        homz.write(sample_id)

        het_cnt=0
        hom_cnt=0
        curr_chr=0
        f = open(fname)
        for l in f.xreadlines():
            lsplit=l.split('\t')
            while lsplit[0] != chromosomes[curr_chr] and curr_chr+1<len(chromosomes):
                hetz.write("\t" + str(het_cnt))
                homz.write("\t" + str(hom_cnt))
                het_cnt=0
                hom_cnt=0
                curr_chr+=1
            # if the call is in a chromosome that is not in the list (or calls could be out of order)
            if lsplit[0] != chromosomes[curr_chr] and curr_chr+1 >= len(chromosomes):
                raise Exception("Unrecognized or out-of-order chromosome: "+lsplit[0])
            
            if lsplit[5] == 'het': het_cnt+=1
            elif lsplit[5] == 'hom': hom_cnt+=1
            else:
                raise Exception("Variant is not a het nor a hom")

        f.close()

        # flush out the counts for the last/remaining chromosomes
        while curr_chr<len(chromosomes):
            hetz.write("\t" + str(het_cnt))
            homz.write("\t" + str(hom_cnt))
            het_cnt=0
            hom_cnt=0             
            curr_chr+=1
        hetz.write('\n')
        homz.write('\n') 

    hetz.close()
    homz.close()
    

@transform(prepare_annovar_inputs, suffix('.avinput'), ['.avinput.variant_function.stats','.avinput.exonic_variant_function.stats'])
def get_stats_on_raw_variants(input, outputs):
    """ annotate functional change in raw variants, get stats, and remove annotated files """
    get_stats_on_prefiltered_variants(input, outputs)


@transform(filter_common_1000genomes, suffix('.hg19_EUR.sites.2012_04_filtered'), 
                                        ['.hg19_EUR.sites.2012_04_filtered.variant_function.stats',
                                         '.hg19_EUR.sites.2012_04_filtered.exonic_variant_function.stats'])
def get_stats_on_1kg_filtered_variants(inputs, outputs):
    """ annotate functional change in 1kg filtered variants, get stats, and remove annotated files """
    get_stats_on_prefiltered_variants(inputs[0], outputs)

# equivalent of annotate_rare_variants (input and output files are the same)
@transform(filter_common_inhouse, suffix('.common_inhouse_filtered'), 
                                        ['.common_inhouse_filtered.variant_function.stats',
                                         '.common_inhouse_filtered.exonic_variant_function.stats'])
def get_stats_on_inhouse_filtered_variants(inputs, outputs):
    """ annotate functional change in inhouse-exomes filtered variants, get stats, and remove annotated files """
    get_stats_on_prefiltered_variants(inputs[0], outputs, cleanup=False)



@merge([get_stats_on_raw_variants, get_stats_on_1kg_filtered_variants, 
        get_stats_on_inhouse_filtered_variants], 'all_samples_exonic_variant_stats.tsv')
def produce_variant_stats_table(infiles, table_file):
    """ produce a table of per-sample counts of different type of exonic variants """

    var_types=['splicing','UTR3','UTR5','intronic','intergenic','exonic']
    
    # split the input files per task
    sample_no = len(infiles)/3
    raw_variant_files = infiles[0:sample_no]
    kg1_filtered_variant_files = infiles[sample_no:sample_no*2]
    inhouse_filtered_variant_files = infiles[sample_no*2:sample_no*3]

    out = open(table_file,'w')
   
    import itertools
    #header = ['sample'] + ['raw_'+t for t in var_types] + ['rare_'+t for t in var_types] + ['raw_synonymous','rare_synonymous']
    filtering_stages = ['raw_','kg1_','inhouse_']
    header = ['sample'] + \
            [f+t for (f,t) in itertools.product(filtering_stages, var_types)] + \
            [s+'_synonymous' for s in filtering_stages]
    out.write(('\t'.join(header))+'\n')
    for i in range(0,sample_no):
        out.write(os.path.basename(raw_variant_files[i][0]).split('.')[0])
        for fname in [raw_variant_files[i][0], kg1_filtered_variant_files[i][0], \
                        inhouse_filtered_variant_files[i][0]]: # exonic variant stats of raw variants and rare variants
            counts = dict.fromkeys(var_types,0)
            f=open(fname)        
            for l in f.xreadlines():
                for var_type in var_types:
                    if l.find(' '+var_type)>0:
                        counts[var_type] += int(l.split()[0])
                        break
            out.write('\t'+'\t'.join([str(counts[t]) for t in var_types]))
            f.close()

        for fname in [raw_variant_files[i][1], kg1_filtered_variant_files[i][1], \
                        inhouse_filtered_variant_files[i][1]]: # synonymous variants stats of raw and rare variants
            f=open(fname)
            found=False
            for l in f.xreadlines():
                if l.find(" synonymous")>0:
                    found=True
                    out.write('\t'+l.split()[0])
                    break     
            if not found: out.write('\t0')       
            f.close()
        out.write('\n')
    out.close()


@follows(count_hetz_and_homz_per_chr, produce_variant_stats_table)
def variant_qc():
    """ empty task aggregating all variant QC tasks """
    pass


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Main logic


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if __name__ == '__main__':
    if options.just_print:
        pipeline_printout(sys.stdout, options.target_tasks, options.forced_tasks,
                            gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                            verbose=options.verbose,
                            checksum_level = 0)

    elif options.flowchart:
        pipeline_printout_graph (   open(options.flowchart, "w"),
                                    # use flowchart file name extension to decide flowchart format
                                    #   e.g. svg, jpg etc.
                                    os.path.splitext(options.flowchart)[1][1:],
                                    options.target_tasks,
                                    options.forced_tasks,
                                        gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                                    no_key_legend   = not options.key_legend_in_graph)
    else:
        pipeline_run(options.target_tasks, options.forced_tasks,
                            multiprocess    = options.jobs,
                            logger          = stderr_logger,
                            verbose         = options.verbose,
                            gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                            checksum_level  = 0)


