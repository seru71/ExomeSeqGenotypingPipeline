#!/usr/bin/env python
"""

    pipeline_multisample.py
                        [--bamdir PATH]
                        [--groups INT]
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


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   options


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


if __name__ == '__main__':
    from optparse import OptionParser
    import StringIO

    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --settings PIPELINE_SETTINGS.CFG --groups NUMBER [more_options]")
    parser.add_option("-s", "--settings", dest="pipeline_settings",
                        metavar="FILE",
                        type="string",
                        help="File containing all the settings for the analysis.")                  
                            
    parser.add_option("-g", "--groups", dest="groups",
                        type="int",
                        default=1,
                        help="Split dataset into smaller groups for multisample snp caling.")



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
    try: input_bams = config.get('Inputs','input-bams')
    except (ConfigParser.NoOptionError):
        sys.stderr.write('No input-bams setting in config file, using fastq-files.')
        input_bams=None
        input_fastqs=config.get('Inputs','input-fastqs')

    # reference files
    reference = config.get('Resources','reference-genome')
    bwa_reference = config.get('Resources','bwa-reference')
    dbsnp = config.get('Resources','dbsnp-vcf')
    hapmap = config.get('Resources','hapmap-vcf')
    omni = config.get('Resources','1000genomes-omni-vcf')
    snps_1kg = config.get('Resources','1000genomes-snps-vcf')
    indels_1kg = config.get('Resources','1000genomes-indels-vcf')
    mills = config.get('Resources','mills-indels-vcf')
    capture = config.get('Resources','capture-regions-bed')
    #capture_qualimap = config.get('Resources','capture-regions-bed-for-qualimap')
    exome = config.get('Resources', 'exome-regions-bed')
    adapters = config.get('Resources', 'trimmomatic-adapters')
    
    # tools 
    trimmomatic = config.get('Tools','trimmomatic-jar')
    bwa = config.get('Tools','bwa-binary')
    samtools = config.get('Tools','samtools-binary')
    java = config.get('Tools','java-binary')
    picard = config.get('Tools','picard-tools-path')
    qualimap = config.get('Tools','qualimap')
    gatk = config.get('Tools','gatk-jar')
    vcftools = config.get('Tools','vcftools')
    # other
    n_cpus = config.getint('Other','n-cpus')





#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   imports


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

from ruffus import *
import subprocess
# import drmaa
import resource

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

def rename(old_file, new_file):
    """rename file"""
    os.rename(old_file, new_file)

def remove(f):
    """remove file"""
    os.remove(f)
                            
def index_bam(bam):
    """Use samtools to create an index for the bam file"""
    run_cmd('%s index %s' % (samtools, bam))
                          
def bam_quality_score_distribution(bam,qs,pdf):
    """Calculates quality score distribution histograms"""
    run_cmd("java -jar {picard}/QualityScoreDistribution.jar \
             CHART_OUTPUT={chart} \
             OUTPUT={output} \
             INPUT={bam} \
             VALIDATION_STRINGENCY=SILENT".format(
                 picard=picard,
                 chart=pdf,
                 output=qs,
                 bam=bam
             ))

def bam_alignment_metrics(bam,metrics):
    """Collects alignment metrics for a bam file"""
    run_cmd("java -jar {picard}/CollectAlignmentSummaryMetrics.jar \
             REFERENCE_SEQUENCE={reference} \
             OUTPUT={output} \
             INPUT={bam} \
             VALIDATION_STRINGENCY=SILENT".format(
                 picard=picard,
                 reference=reference,
                 output=metrics,
                 bam=bam
             ))
    
def bam_coverage_metrics(input_bam, output):
    """ Calculates and outputs bam coverage statistics """
    run_cmd("{java} -Xmx4g -jar {gatk} \
            -R {reference} \
            -T DepthOfCoverage \
            -o {output} \
            -I {input} \
            -L {capture} \
            -ct 8 -ct 20 -ct 30 \
            --omitDepthOutputAtEachBase --omitLocusTable \
            ".format(java=java, gatk=gatk,
                reference=reference,
                output=output,
                input=input_bam,
                capture=capture
            ))
 
# def qualimap_bam(input_bam, output_dir):
#     """ Generates Qualimap bam QC report """
#     # create necessary folders first
#     if not os.path.exists('qc'): os.mkdir('qc')
#     if not os.path.exists('qc/qualimap/'): os.mkdir('qc/qualimap')
#     if not os.path.exists(output_dir): os.mkdir(output_dir)
#     run_cmd("{qualimap} bamqc -bam {bam} \
#              -c -outformat PDF \
#              -gff {target} \
#              -gd HUMAN -os \
#              -outdir {dir} &> {dir}/qualimap.err".format(
#                 qualimap=qualimap,
#                 bam=input_bam,
#                 target=capture_qualimap,
#                 dir=output_dir))

    
    
#
#
# The two functions below seem redundant
#
#    
#            
#def merge_indel_and_snp_vcf(snp,indel, output):
#    """Merges vcf files from the batch run"""
#    run_cmd("{java} -Xmx4g -jar {gatk} \
#            -R {reference} \
#            -T CombineVariants \
#            -V:SNP {snp} \
#            -V:INDEL {indel} \
#            -o {output}".format(
#                java=java,
#                gatk=gatk,
#                reference=reference,
#                snp=snp,
#                indel=indel,
#                output=output
#            )) 
#
#def variant_annotator(bam, vcf, output):
#    run_cmd("{java} -Xmx16g -jar {gatk} \
#            -R {reference} \
#            -T VariantAnnotator \
#            -I {bam} \
#            -o {output} \
#            -A Coverage \
#            --variant {vcf} \
#            --dbsnp {dbsnp} \
#            ".format(java=java, gatk=gatk,
#                bam=bam,
#                reference=reference,
#                output=output,
#                vcf=vcf,
#                dbsnp=dbsnp
#            )) 


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

#   Pipeline


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#       Put pipeline code here

def get_sample_ids():
    if input_bams != None:
        files = glob.glob(input_bams)
        return [ os.path.splitext(os.path.basename(f))[0] for f in files ]
    else:
        fqs = sorted(glob.glob(input_fastqs))
        return [ os.path.basename(fqs[i]).split('_L')[0] for i in range(0,len(fqs),2) ]

def get_num_files():
    if input_bams != None: 
        return len(glob.glob(input_bams))
    else: 
        return len(glob.glob(input_fastqs)) / 2


def generate_parameters():
    fqs = sorted(glob.glob(input_fastqs))
    parameters = []
    for i in range(0,len(fqs),2):
        fq1=fqs[i]
        fq2=fqs[i+1]
        prefix = os.path.basename(fq1).split('_L')[0]
        parameters.append([[fq1, fq2], 
                          [os.path.join(prefix,prefix+'.clean.fq1.gz'), os.path.join(prefix,prefix+'.clean.fq2.gz')],
                          prefix])
    for job_parameters in parameters:
            yield job_parameters



#
#
# Prepare input directory structure
#

@files(generate_parameters)
def trim_reads(fastqs_in, fastqs_out, dirname):
    """ Create sample dir, and trim the reads """
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    
    run_cmd("{java} -jar {trm} PE -phred33 {fq1} {fq2} {fq1_clean} {fq1_drop} {fq2_clean} {fq2_drop} \
            ILLUMINACLIP:{adapters}:2:30:10 \
            LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50".format(
                java=java,
                trm=trimmomatic,
                fq1=fastqs_in[0], fq2=fastqs_in[1],
                fq1_clean=fastqs_out[0], fq2_clean=fastqs_out[1],
                fq1_drop=fastqs_out[0]+".drop",
                fq2_drop=fastqs_out[1]+".drop",
                adapters=adapters))

    # rm dropped reads
    remove(fastqs_out[0]+".drop")
    remove(fastqs_out[1]+".drop")


#@transform(trim_reads, regex(".clean.fq[12].gz"), ".bam")
@transform(trim_reads, formatter(".*/(?P<SAMPLE_ID>[^/]+).clean.fq1.gz", None), "{SAMPLE_ID[0]}/{SAMPLE_ID[0]}.bam","{SAMPLE_ID[0]}")
def align(fastqs, bam, sample_name):
    """ Align the reads to the reference and sort """
    
    # align
    run_cmd("{bwa} mem {ref} {fq1} {fq2} | {samtools} view -bT {ref} - | {samtools} sort -f - {bam}".format(
            bwa=bwa,
            samtools=samtools,
            ref=bwa_reference,
            fq1=fastqs[0],
            fq2=fastqs[1],
            bam=bam))
 
   
    # add read groups...
    run_cmd("{java} -jar {add_rg} I={bam} O={rg_bam} LB={sample} PL=illumina PU=unit SM={sample}".format( 
             java=java, 
             add_rg=os.path.join(picard,"AddOrReplaceReadGroups.jar"),
             bam=bam,
             rg_bam=bam+".tmp",
             sample=sample_name))
    # ..."in-place"
    rename(bam+".tmp", bam)
    

@transform(align, suffix(".bam"), '.bam.bai')
def index(bam, output):
    """create bam index"""
    index_bam(bam)


#
#
# QC the raw bam files
#

@follows(index)
@transform(align, suffix(".bam"), '.quality_score')
def qc_raw_bam_quality_score_distribution(input_bam, output):
    """docstring for metrics1"""
    bam_quality_score_distribution(input_bam, output, output + '.pdf')

@follows(index)
@transform(align, suffix(".bam"), '.metrics')
def qc_raw_bam_alignment_metrics(input_bam, output):
    """docstring for metrics1"""
    bam_alignment_metrics(input_bam, output)
    
@follows(index)
@transform(align, suffix(".bam"), '.coverage.sample_summary', r'\1.coverage')
def qc_raw_bam_coverage_metrics(input_bam, output, output_format):
    bam_coverage_metrics(input_bam, output_format)

@follows(index)
@transform(align, formatter(".*/(?P<SAMPLE_ID>[^/]+).bam"), '{subpath[0][1]}/qc/qualimap/{SAMPLE_ID[0]}')
def qc_raw_bam_qualimap_report(input_bam, output_dir):
    qualimap_bam(input_bam, output_dir)

#@follows(qc_raw_bam_quality_score_distribution, qc_raw_bam_alignment_metrics, qc_raw_bam_coverage_metrics)
@follows(qc_raw_bam_coverage_metrics, qc_raw_bam_qualimap_report)
def raw_bam_qc():
    """ Aggregates raw bam quality control steps """
    pass



    

#
#
# Prepare the bam files for variant calling
#

#
# don't remove duplicates - this is amplicon data
# ============================================================
#
#def dup_removal_picard(bam,output):
#    """Use Picard to remove duplicates"""
#    run_cmd('%s -Xmx4096m -jar %s/CleanSam.jar \
#             INPUT=%s \
#              OUTPUT=%s \
#              VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR CREATE_INDEX=TRUE' 
#              % (java, picard, bam, output))
#    
# def dup_mark_picard(bam,output):
#     """Use Picard to mark duplicates"""
#     run_cmd('%s -Xmx4096m -jar %s/MarkDuplicates.jar \
#             TMP_DIR=/export/astrakanfs/stefanj/tmp \
#             REMOVE_DUPLICATES=true \
#             INPUT=%s \
#             OUTPUT=%s \
#             METRICS_FILE=%s.dup_metrics \
#             VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR CREATE_INDEX=true'
#             % (java, picard, bam, output, bam))
# 
# 
# @follows(index)
# @transform(link, suffix(".bam"), '.dedup.bam')
# def remove_dups(bam, output):
#     """Use samtools for dupremoval"""
#     #run_cmd('samtools rmdup %s %s' % (bam,output))
#     dup_mark_picard(bam, output)
#     # remove(input)
# 
# 
# #@follows(remove_dups)
# @transform(remove_dups, suffix(".dedup.bam"), '.dedup.bam.bai')
# def index_dups(bam, output):
#     """create bam index"""
#     index_bam(bam)


@follows(index)
@transform(align, suffix(".bam"), '.realign.intervals')
def find_realignment_intervals(input_bam, intervals):
    """Find regions to be re-aligned due to indels"""
    run_cmd("%s -Xmx4g -jar %s \
             -T RealignerTargetCreator \
             -I %s \
             -R %s \
             -L %s \
             -known %s \
             -known %s \
             -o %s \
             -dcov 3000 -nt %d" 
             % (java, gatk, input_bam, reference, exome, indels_1kg, mills, intervals, n_cpus))


#@follows(find_realignment_intervals)
@transform(find_realignment_intervals, suffix(".realign.intervals"), '.realigned.bam', r'\1.bam')
def indel_realigner(intervals_file, realigned_bam, input_bam):
    """Re-aligns regions around indels"""
    run_cmd("%s -Xmx4g -jar %s \
             -T IndelRealigner \
             -I %s \
             -R %s \
             -L %s \
             -targetIntervals %s \
             -known %s \
             -known %s \
             -o %s" 
             % (java, gatk, input_bam, reference, exome, intervals_file, indels_1kg, mills, realigned_bam))


#@follows(indel_realigner)
@transform(indel_realigner, suffix('.realigned.bam'), '.gatk.bam.recal_data.grp')
def recalibrate_baseq1(input_bam, output):
    """Base quality score recalibration in bam file - first pass """
    index_bam(input_bam)
    run_cmd("%s -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx6g -jar %s \
            -T BaseRecalibrator \
            -R %s \
            -L %s \
            -knownSites %s \
            -I %s \
            -o %s \
            -nct %d"
            % (java, gatk, reference, exome, dbsnp, input_bam, output, n_cpus))
    

# This custom check ensures that the recalibrate_baseq2 step is not run in --rebuild_mode if the .gatk.bam exists
# This way metric_coverage target can be run after the intermediate files (e.g. .realigned.bam) are removed
def is_gatk_bam_missing(inputs, gatk_bam):
    # ignore input files of this step: input[0] - .realigned.bam, and input[1] - .recal_data.grp file
    if not os.path.exists(gatk_bam):
        return True, "Missing file %s" % gatk_bam
    else:
        return False, "File %s exists" % gatk_bam

@follows(recalibrate_baseq1)
@transform(indel_realigner, suffix('.realigned.bam'), add_inputs(r'\1.gatk.bam.recal_data.grp'),'.gatk.bam')
@check_if_uptodate(is_gatk_bam_missing)
def recalibrate_baseq2(inputs, output_bam):
    """Base quality score recalibration in bam file
        Part 2: rewrite quality scores into a new bam file"""   
    bam = inputs[0]
    recal_data = inputs[1]
    run_cmd("%s -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx4g -jar %s \
            -T PrintReads \
            -R %s \
            -L %s \
            -I %s \
            --out %s \
            -BQSR %s \
            -nct %d" 
            % (java, gatk, reference, exome, bam, output_bam, recal_data, n_cpus) )
  
    # remove(inputs[0])


#
#
# gatk.bam-level QC measurements
#

#@follows(recalibrate_baseq2)
@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.quality_score')
def qc_gatk_bam_quality_score_distribution(input_bam, output):
    """docstring for metrics1"""
    bam_quality_score_distribution(input_bam, output, output + '.pdf')

#@follows(recalibrate_baseq2)
@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.metrics')
def qc_gatk_bam_alignment_metrics(input_bam, output):
    """docstring for metrics1"""
    bam_alignment_metrics(input_bam, output)

#@follows(recalibrate_baseq2)
@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.coverage.sample_summary', r'\1.coverage')
def qc_gatk_bam_coverage_metrics(input_bam, output, output_format):
    bam_coverage_metrics(input_bam, output_format)

@transform(recalibrate_baseq2, formatter(".*/(?P<SAMPLE_ID>[^/]+).bam"), '{path[0]}/qc/qualimap/{SAMPLE_ID[0]}')
def qc_gatk_bam_qualimap_report(input_bam, output_dir):
    qualimap_bam(input_bam, output_dir)


#@follows(qc_gatk_bam_quality_score_distribution, qc_gatk_bam_alignment_metrics, qc_gatk_bam_coverage_metrics, qc_gatk_bam_qualimap_report)
@follows(qc_gatk_bam_coverage_metrics, qc_gatk_bam_qualimap_report)
def gatk_bam_qc():
    """ Aggregates gatk_bam quality control steps """
    pass



#
#
# Bam reduce & (batch) variant calling using UnifiedGenotyper
#


@jobs_limit(6)
#@follows(recalibrate_baseq2)
@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.reduced.bam')
def reduce_bam(bam, output):
    """Reduces the BAM file using read based compression that keeps only essential information for variant calling"""
    run_cmd("%s -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx6g -jar %s \
            -T ReduceReads \
            -R %s \
            -I %s \
            -o %s"
            % (java, gatk, reference, bam, output))


def split_seq(seq, num_pieces):
    """ split a list into pieces passed as param """
    start = 0
    for i in xrange(num_pieces):
        stop = start + len(seq[i::num_pieces])
        yield seq[start:stop]
        start = stop

def multisample_variant_call(bams, output):
    """Perform multi-sample variant calling using GATK"""
    cmd = "nice %s -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx8g -jar %s \
            -T UnifiedGenotyper \
            -R %s \
            -o %s \
            -glm BOTH \
            -nt %s \
            --dbsnp %s " % (java, gatk, reference, output, options.jobs, dbsnp)
    for bam in bams:
        cmd = cmd + '-I {} '.format(bam)
    #log the results
    cmd = cmd + '&> {}.log'.format(output)
    run_cmd(cmd)
    
def merge_batch_vcf(vcfs, output):
    """Merges vcf files from the batch run"""
    if len(vcfs) == 1:
        run_cmd('cp {vcf} {output}'.format(vcf = vcfs[0], output = output))
    else:
        merging = ''
        for i in range(len(vcfs)):
            merging = merging + ' -V:batch{number} {file}'.format(number=i,file=vcfs[i])
        run_cmd("{java} -Xmx4g -jar {gatk} \
                -R {reference} \
                -T CombineVariants \
                -o {output} \
                {files}".format(
                    java=java,
                    gatk=gatk,
                    reference=reference,
                    output=output,
                    files=merging
                )) 

@merge(reduce_bam, 'multisample.gatk.vcf')
def call_variants(infiles, output):
    """Splits the files into s number of batches, calls the variants, and merges them back"""
    #splits the calculations into batches
    counter = 1
    batches = []
    for batch in split_seq(infiles, options.groups):
        multisample_variant_call(batch, 'batch{}.vcf'.format(counter))
        batches.append('batch{}.vcf'.format(counter))
        counter += 1
    merge_batch_vcf(batches,output)
    for batch in batches:
        remove(batch)
        remove(batch + '.idx')


# @follows('call_variants')
# @files('multisample.gatk.vcf', ['multisample.gatk.snp.vcf','multisample.gatk.indel.vcf'])
# def split_snps_and_indels(vcf, outputs):
#     """Separates the vcf file into only snp calls"""
#    run_cmd("{} -Xmx2g -jar {} -R {} -T SelectVariants \
#            --variant {} -selectType SNP\
#            -o {}".format(java, gatk,reference,vcf, outputs[0]))
#    run_cmd("{} -Xmx2g -jar {} -R {} -T SelectVariants \
#            --variant {} -selectType INDEL\
#            -o {}".format(java, gatk,reference,vcf,outputs[1]))


#
#
# Genotyping using HaplotypeCaller
#

@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.gvcf')
def call_haplotypes(bam, output_gvcf):
    """Perform variant calling using GATK HaplotypeCaller"""
    cmd = "nice %s -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx6g -jar %s \
            -T HaplotypeCaller \
            -R %s \
            -I %s \
            -o %s \
            --emitRefConfidence GVCF \
            --variant_index_type LINEAR \
            --variant_index_parameter 128000 \
            -minPruning 4 \
            -L %s \
            -dcov 3000 \
            --dbsnp %s " % (java, gatk, reference, bam, output_gvcf, exome, dbsnp)
#             -nct %s \
#             -stand_call_conf 50.0 \

    #log the results
    #cmd = cmd + '&> {}.log'.format(output_gvcf)
    run_cmd(cmd)


@merge(call_haplotypes, 'multisample.gatk.vcf')
def genotype_gvcfs(gvcfs, output):
    """Combine the per-sample GVCF files and genotype""" 
    cmd = "nice %s -Xmx16g -jar %s -T GenotypeGVCFs \
            -R %s -o %s -dcov 3000 -nt %s" % (java, gatk, reference, output, options.jobs)
       
    for gvcf in gvcfs:
        cmd = cmd + " --variant {}".format(gvcf)
    
    cmd = cmd + '&> {}.log'.format(output)
    run_cmd(cmd)



#
#
# Variant recalisbration and filtering 
# (not applicable for such a small target)
# 
#===============================================================================
# #@follows(call_variants)
# @follows(genotype_gvcfs)
# @files('multisample.gatk.vcf', ['multisample.gatk.snp.model','multisample.gatk.snp.model.tranches'])
# def find_snp_tranches_for_recalibration(vcf,outputs):
#     """Runs VariantRecalibrator on snps"""
#     cmd = "{java} -Xmx16g -jar {gatk} \
#             -T VariantRecalibrator \
#             -R {reference}  \
#             -input {input} \
#             -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} \
#             -resource:omni,known=false,training=true,truth=true,prior=12.0 {omni} \
#             -resource:1000G,known=false,training=true,truth=false,prior=10.0 {snps_1kg} \
#             -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} \
#             -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an DP \
#             -mode SNP \
#             -recalFile {output} \
#             -tranchesFile {tranches} \
#             -rscriptFile {plots}\
#             -nt {num_jobs}".format(
#                 java=java,
#                 gatk=gatk,
#                 reference=reference,
#                 hapmap=hapmap,
#                 omni=omni,
#                 snps_1kg=snps_1kg,
#                 dbsnp=dbsnp,
#                 input=vcf,
#                 output=outputs[0],
#                 tranches=outputs[1],
#                 plots=outputs[0]+'.plots.R',
#                 num_jobs=options.jobs
#             )
#     if get_num_files() > 10:
#         cmd += " -an InbreedingCoeff"
#     run_cmd(cmd)
# 
# 
# @follows(find_snp_tranches_for_recalibration)
# @files('multisample.gatk.vcf', ['multisample.gatk.indel.model','multisample.gatk.indel.model.tranches'])
# def find_indel_tranches_for_recalibration(vcf,outputs):
#     """Runs VariantRecalibrator on indels"""
#     cmd = "{java} -Xmx16g -jar {gatk} \
#             -T VariantRecalibrator \
#             -R {reference}  \
#             -input {input} \
#             --maxGaussians 4 \
#             -resource:mills,known=false,training=true,truth=true,prior=12.0 {mills} \
#             -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} \
#             -an QD -an DP -an FS -an ReadPosRankSum -an MQRankSum \
#             -mode INDEL \
#             -recalFile {output} \
#             -tranchesFile {tranches} \
#             -rscriptFile {plots}\
#             -nt {num_jobs}".format(
#                 java=java,
#                 gatk=gatk,
#                 reference=reference,
#                 mills=mills,
#                 dbsnp=dbsnp,
#                 input=vcf,
#                 output=outputs[0],
#                 tranches=outputs[1],
#                 plots=outputs[0]+'.plots.R',
#                 num_jobs=options.jobs
#             )
#     if get_num_files() > 10:
#         cmd += " -an InbreedingCoeff"
#     run_cmd(cmd)
# 
# 
# def apply_recalibration_to_snps_or_indels(vcf,recal,tranches,output,mode='SNP'):
#     """Apply the recalibration tranch file to either snps or indels (depending on mode=SNP or mode=INDEL"""
#     run_cmd("{java} -Xmx16g -jar {gatk} \
#             -T ApplyRecalibration \
#             -R {reference} \
#             -input {vcf} \
#             --ts_filter_level 99.9 \
#             -tranchesFile {tranches}  \
#             -recalFile {recal} \
#             -mode {mode} \
#             -nt {num_jobs} \
#             -o {output}".format(
#                 java=java,
#                 gatk=gatk,
#                 reference=reference,
#                 vcf=vcf,
#                 tranches=tranches,
#                 recal=recal,
#                 mode=mode,
#                 num_jobs=options.jobs,
#                 output=output
#             ))
# 
# 
# @follows(find_indel_tranches_for_recalibration)
# @files(['multisample.gatk.vcf','multisample.gatk.snp.model','multisample.gatk.snp.model.tranches'],'multisample.gatk.recalibratedSNPS.rawIndels.vcf')
# def apply_recalibration_filter_snps(inputs,output):
#     apply_recalibration_to_snps_or_indels(inputs[0],inputs[1],inputs[2],output,mode='SNP')
#     # remove(input[0])
#     #     remove(input[1])
#     #     remove(input[2])
# 
# 
# @follows(apply_recalibration_filter_snps)
# @files(['multisample.gatk.recalibratedSNPS.rawIndels.vcf','multisample.gatk.indel.model','multisample.gatk.indel.model.tranches'],'multisample.gatk.preHardFiltering.vcf')
# def apply_recalibration_filter_indels(inputs,output):
#     apply_recalibration_to_snps_or_indels(inputs[0],inputs[1],inputs[2],output,mode='INDEL')
#===============================================================================

#@follows(apply_recalibration_filter_indels)
#@files('multisample.gatk.preHardFiltering.vcf', 'multisample.gatk.markedHardFiltering.vcf')
@transform(genotype_gvcfs, suffix('.vcf'), '.markedHardFiltering.vcf')
def filter_variants(input_vcf, output_vcf):
    """docstring for apply_indel_filter"""
    run_cmd('{java} -Xmx12g -jar {gatk} \
            -T VariantFiltration \
            -o {output} \
            --variant {input} \
            --filterExpression "QD < 3.0" \
            --filterExpression "DP < 6" \
            --filterName QDFilter   \
            --filterName DPFilter   \
            -R {reference}'.format(
                java=java,
                gatk=gatk,
                output=output_vcf,
                input=input_vcf,
                reference=reference
            ))
    # remove(input)


#@follows(filter_variants)
#@files('multisample.gatk.markedHardFiltering.vcf', 'multisample.gatk.analysisReady.vcf')
@transform(filter_variants, suffix('.markedHardFiltering.vcf'), '.analysisReady.vcf')
def remove_filtered(input_vcf, output_vcf):
    """Remove filtered variants"""
    run_cmd("{java} -Xmx12g -jar {gatk} \
            -T SelectVariants \
            -R {reference} \
            --variant {input} \
            -o {output} \
            -env -ef".format(
                java=java,
                gatk=gatk,
                reference=reference,
                input=input_vcf,
                output=output_vcf
            ))
    # remove(input)


#@follows(remove_filtered)
@transform(remove_filtered, suffix('.analysisReady.vcf'),'.analysisReady.exome.vcf')
def final_calls(vcf, output):
    """ Produce the final variant calls in the exome regions """
    output = output[:-10]
    # apply filters to the vcf file to limit calling to exome region    
    run_cmd("%s --vcf %s \
             --out %s \
             --recode \
             --bed %s \
             --keep-INFO-all"
            % (vcftools, vcf, output, exome))
    rename('multisample.gatk.analysisReady.recode.vcf','multisample.gatk.analysisReady.exome.vcf')


def split_snp_parameters():
    exome_vcf = 'multisample.gatk.analysisReady.exome.vcf'
    for s_id in get_sample_ids():
        yield [exome_vcf, s_id + '/' + s_id + '.exome.vcf', s_id]

def cleanup_files():
    run_cmd("rm -rf */*.recal_data.csv */*.realign* */*.dedup* */*.log *.log \
            *.to_filter* multisample.gatk.snp.recal batch* \
            multisample.gatk.recalibratedSNPS.rawIndels.vcf \
            multisample.gatk.indel.model.* multisample.gatk.snp.model.* \
            multisample.gatk.analysisReady.vcf.vcfidx \
            multisample.gatk.analysisReady.vcf.idx \
            multisample.gatk.recalibratedSNPS.rawIndels.vcf.idx \
            ")
#            multisample.gatk.analysisReady.vcf \


@posttask(cleanup_files)
@follows('final_calls')
@files(split_snp_parameters)
def split_snps(vcf, output, sample):
    """ Split variants by sample, and use sample-specific statistics to filter: AD and DP"""
    AD_threshold=5
    DP_threshold=8
    # what if there is multiple alt alleles??? is the first one most covered      
    run_cmd("{java} -Xmx2g -jar {gatk} -R {ref} \
            -T SelectVariants \
            --variant {vcf} \
            -sn {sample} \
            -select 'vc.getGenotype(\"{sample}\").getAD().1 >= {ad_thr} && vc.getGenotype(\"{sample}\").getDP() >= {dp_thr}' \
            -o {out} \
            ".format(java=java, gatk=gatk, ref=reference, 
                vcf=vcf, sample=sample, out=output, 
                ad_thr=AD_threshold, 
                dp_thr=DP_threshold))


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

