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
import string

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   user definable options


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#path to binaries
script_path = os.path.dirname(os.path.relpath(__file__))
java = os.path.join(script_path,'../src/jre1.7.0/bin/java')
picard = os.path.join(script_path,'../src/picard-tools')
gatk = os.path.join(script_path,'../src/GenomeAnalysisTK/GenomeAnalysisTK.jar')
#reference files
reference = os.path.join(script_path,'../reference/human_g1k_v37.clean.fasta')
dbsnp = os.path.join(script_path,'../reference/dbsnp_137.b37.vcf')
hapmap = os.path.join(script_path,'../reference/hapmap_3.3.b37.sites.vcf')
omni = os.path.join(script_path,'../reference/1000G_omni2.5.b37.sites.vcf')
mills = os.path.join(script_path,'../reference/Mills_and_1000G_gold_standard.indels.b37.vcf')
#capture = os.path.join(script_path,'../reference/Nimblegen_SeqCap_EZ_Exome_v2_37_targetRegOnly_g1k.bed')
#exome = os.path.join(script_path,'../reference/Nimblegen_SeqCap_EZ_Exome_v2_37_targetRegOnly_wingspan_g1k.bed')

n_cpus = 2


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   options


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


if __name__ == '__main__':
    from optparse import OptionParser
    import StringIO

    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --bamdir BAM_DIR --target_regions EXOME.BED --groups NUMBER [more_options]")
    parser.add_option("-b", "--bamdir", dest="bam_dir",
                        metavar="FILE",
                        type="string",
                        help="Directory containing all the bams for analysis.")
                        
#    parser.add_option("-r","--reference",dest="reference",
#                      metavar="FILE", type="string",
#                      default="../reference/human_g1k_v37.clean.fasta",
#                      help="Fasta file with the reference genome")                 
                     
    parser.add_option("--cr","--capture_regions", dest="capture_regions",
                      metavar="FILE", 
                      type="string",
                      default="../reference/Nimblegen_SeqCap_EZ_Exome_v2_37_targetRegOnly_g1k.bed", 
                      help="Bed file with capture regions (used for depth calculations)")                    

    parser.add_option("--er","--exome_regions", dest="exome_regions",
                      metavar="FILE", 
                      type="string", 
                      default="../reference/Nimblegen_SeqCap_EZ_Exome_v2_37_targetRegOnly_wingspan_g1k.bed",
                      help="Bed file with exome regions (used for filtering; can be wider than capture regions)")                    
                            
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
    mandatory_options = ['bam_dir']

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


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   imports


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

from ruffus import *
import subprocess
# import drmaa
import resource

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Functions and variables


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

def rename(file,new_file):
    """rename file"""
    os.rename(file,new_file)

def remove(file):
    """remove file"""
    os.remove(file)
                            
def split_seq(seq, num_pieces):
    """ split a list into pieces passed as param """
    start = 0
    for i in xrange(num_pieces):
        stop = start + len(seq[i::num_pieces])
        yield seq[start:stop]
        start = stop

def index_bam(bam):
    """Use samtools to create an index for the bam file"""
    run_cmd('samtools index %s' % bam)
    
def find_realigns(bam, intervals_file):
    """Finds regions needing realignment"""
    run_cmd("%s -Xmx4g -jar %s \
             -T RealignerTargetCreator \
             -I %s \
             -R %s \
             -o %s \
             -nt %d" 
             % (java, gatk, bam, reference, intervals_file, n_cpus))
             
def run_realigner(bam, intervals_file, realigned_bam):
    """Performs local realignment of reads based on misalignments due to the presence of indels"""
    run_cmd("%s -Xmx4g -jar %s \
             -T IndelRealigner \
             -I %s \
             -R %s \
             -targetIntervals %s \
             -o %s" 
             % (java, gatk, bam, reference, intervals_file, realigned_bam))

# def dup_removal_picard(bam,output):
#     """Use Picard to remove duplicates"""
#     run_cmd('java -Xmx4096m -jar %s/CleanSam.jar INPUT=%s OUTPUT=%s VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR CREATE_INDEX=TRUE' 
#                 % (picard,bam,output))

 
def dup_removal_samtools(bam,output):
    """Use samtools for dupremoval"""
    run_cmd('samtools rmdup %s %s' % (bam,output))
    
# def dup_mark_picard(bam,output):
#     """Use Picard to mark duplicates"""
#     run_cmd('java -Xmx4096m -jar %s/MarkDuplicates.jar TMP_DIR=/export/astrakanfs/stefanj/tmp REMOVE_DUPLICATES=true INPUT=%s OUTPUT=%s METRICS_FILE=%s.dup_metrics VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR CREATE_INDEX=true'
#                 % (picard,bam,output,bam))

def base_recalibrator(bam, recal_data):
    """First pass of the recalibration step"""
    run_cmd("%s -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx8g -jar %s \
            -T BaseRecalibrator \
            -R %s \
            -knownSites %s \
            -I %s \
            -o %s \
            -nct %d"
            % (java, gatk, reference, dbsnp, bam, recal_data, n_cpus))

def print_recalibrated(bam, recal_data, output):
    """uses GATK to rewrite quality scores using the recal_data"""
    run_cmd("%s -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx4g -jar %s \
            -T PrintReads \
            -R %s \
            -I %s \
            --out %s \
            -BQSR %s \
            -nct %d" 
            % (java, gatk, reference, bam, output, recal_data, n_cpus) )
            
# def count_covariates(bam, recal_data):
#     """Uses GATK to count covariates"""
#     run_cmd("java -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx4g -jar %s \
#             -T CountCovariates \
#             -l INFO \
#             -R %s \
#             -knownSites %s \
#             --default_platform illumina \
#             -cov ReadGroupCovariate \
#             -cov QualityScoreCovariate \
#             -cov CycleCovariate \
#             -cov DinucCovariate \
#             -I %s \
#             -recalFile %s"
#             % (gatk, reference, dbsnp, bam, recal_data))
            
# def table_recalibration(bam, recal_data, output):
#     """uses GATK to rewrite quality scores using the recal_data"""
#     run_cmd("java -Xmx4g -jar %s \
#             -T TableRecalibration \
#             --default_platform illumina \
#             -R %s \
#             --preserve_qscores_less_than 5 \
#             -l INFO \
#             -I %s \
#             --out %s \
#             -recalFile %s" 
#             % (gatk, reference, bam, output, recal_data) )

def reduce_reads(bam, output):
    """Reduces the BAM file using read based compression that keeps only essential information for variant calling"""
    run_cmd("%s -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx15g -jar %s \
            -T ReduceReads \
            -R %s \
            -I %s \
            -o %s"
            % (java, gatk, reference, bam, output))
            

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

def filter_by_exome_region(vcf, output):
    """Apply filters to the vcf file to limit calling to exome region"""
    
    exome = options.exome_regions 
    if exome==None or exome=="": exome = options.capture_regions 
    
    run_cmd("vcftools --vcf %s \
             --out %s \
             --recode \
             --bed %s \
             --keep-INFO-all "
            % (vcf, output, exome))

def multisample_variant_call(files, output):
    """Perform multi-sample variant calling using GATK"""
    cmd = "nice %s -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx8g -jar %s \
            -T UnifiedGenotyper \
            -R %s \
            -o %s \
            -glm BOTH \
            -nt %s \
            --dbsnp %s " % (java, gatk, reference, output, options.jobs, dbsnp)
    for file in files:
        cmd = cmd + '-I {} '.format(file)
    #log the results
    cmd = cmd + '&> {}.log'.format(output)
    run_cmd(cmd)

def multisample_haplotype_call(files, output):
    """Perform multi-sample variant calling using GATK"""
    cmd = "nice %s -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx24g -jar %s \
            -T HaplotypeCaller \
            -R %s \
            -o %s \
            -minPruning 4 \
            -stand_call_conf 50.0 \
            -stand_emit_conf 30.0 \
            -L %s \
            -nct %s \
            --dbsnp %s " % (java, gatk, reference, output, capture, options.jobs, dbsnp)
    for file in files:
        cmd = cmd + '-I {} '.format(file)
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
            
def merge_indel_and_snp_vcf(snp,indel, output):
    """Merges vcf files from the batch run"""
    run_cmd("{java} -Xmx4g -jar {gatk} \
            -R {reference} \
            -T CombineVariants \
            -V:SNP {snp} \
            -V:INDEL {indel} \
            -o {output}".format(
                java=java,
                gatk=gatk,
                reference=reference,
                snp=snp,
                indel=indel,
                output=output
            )) 

#
# what if there is multiple alt alleles??? is the first one most covered  
#
def split_variants_by_sample(vcf, sample, output, ad_threshold=5, dp_threshold=8):
    
    run_cmd("{java} -Xmx2g -jar {gatk} -R {ref} \
            -T SelectVariants \
            --variant {vcf} \
            -sn {sample} \
            -select 'vc.getGenotype(\"{sample}\").getAD().1 >= {ad_thr} && vc.getGenotype(\"{sample}\").getDP() >= {dp_thr}' \
            -o {out} \
            ".format(java=java, gatk=gatk, ref=reference, 
                vcf=vcf, sample=sample, out=output, 
                ad_thr=ad_threshold, 
                dp_thr=dp_threshold))

def split_snps(vcf,snp_file):
    """Select for snps in the vcf file"""
    run_cmd("{} -Xmx2g -jar {} -R {} -T SelectVariants \
            --variant {} -selectType SNP\
            -o {}".format(java, gatk,reference,vcf,snp_file))
            
def split_indels(vcf,indel_file):
    """Select for indels in the vcf file"""
    run_cmd("{} -Xmx2g -jar {} -R {} -T SelectVariants \
            --variant {} -selectType INDEL\
            -o {}".format(java, gatk,reference,vcf,indel_file))

def variant_annotator(bam, vcf, output):
    run_cmd("{java} -Xmx16g -jar {gatk} \
            -R {reference} \
            -T VariantAnnotator \
            -I {bam} \
            -o {output} \
            -A Coverage \
            --variant {vcf} \
            --dbsnp {dbsnp} \
            ".format(java=java, gatk=gatk,
                bam=bam,
                reference=reference,
                output=output,
                vcf=vcf,
                dbsnp=dbsnp
            ))

def recalibrate_snps(vcf,output):
    """Runs VariantRecalibrator on the snps file"""
    cmd = "{java} -Xmx16g -jar {gatk} \
            -T VariantRecalibrator \
            -R {reference}  \
            -input {input} \
            --maxGaussians 6 \
            -percentBad 0.01 -minNumBad 1000 \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} \
            -resource:omni,known=false,training=true,truth=false,prior=12.0 {omni} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} \
            -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
            -mode SNP \
            -recalFile {output} \
            -tranchesFile {tranches} \
            -rscriptFile {plots}\
            -nt {num_jobs}".format(
                java=java,
                gatk=gatk,
                reference=reference,
                hapmap=hapmap,
                omni=omni,
                dbsnp=dbsnp,
                input=vcf,
                output=output,
                tranches=output+'.tranches',
                plots=output+'.plots.R',
                num_jobs=options.jobs
            )
    if get_num_files() > 10:
        cmd += " -an InbreedingCoeff"
    run_cmd(cmd)

def recalibrate_indels(vcf,output):
    """Runs VariantRecalibrator on the indels file"""
    cmd = "{java} -Xmx16g -jar {gatk} \
            -T VariantRecalibrator \
            -R {reference}  \
            -input {input} \
            --maxGaussians 4 \
            -percentBad 0.01 \
            -minNumBad 1000 \
            -resource:mills,known=false,training=true,truth=true,prior=12.0 {mills} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} \
            -an DP -an FS -an ReadPosRankSum -an MQRankSum \
            -mode INDEL \
            -recalFile {output} \
            -tranchesFile {tranches} \
            -rscriptFile {plots}\
            -nt {num_jobs}".format(
                java=java,
                gatk=gatk,
                reference=reference,
                mills=mills,
                dbsnp=dbsnp,
                input=vcf,
                output=output,
                tranches=output+'.tranches',
                plots=output+'.plots.R',
                num_jobs=options.jobs
            )
    if get_num_files() > 10:
        cmd += " -an InbreedingCoeff"
    run_cmd(cmd)

def apply_recalibration_snps(vcf,recal,tranches,output):
    """Apply the recalibration tranch file calculated in recalibrate_snps"""
    run_cmd("{java} -Xmx16g -jar {gatk} \
            -T ApplyRecalibration \
            -R {reference} \
            -input {vcf} \
            --ts_filter_level 99.9 \
            -tranchesFile {tranches}  \
            -recalFile {recal} \
            -mode SNP \
            -nt {num_jobs} \
            -o {output}".format(
                java=java,
                gatk=gatk,
                reference=reference,
                vcf=vcf,
                tranches=tranches,
                recal=recal,
                output=output,
                num_jobs=options.jobs
            ))

def apply_recalibration_indels(vcf,recal,tranches,output):
    """Apply the recalibration tranch file calculated in recalibrate_snps"""
    run_cmd("{java} -Xmx16g -jar {gatk} \
            -T ApplyRecalibration \
            -R {reference} \
            -input:name,VCF {vcf} \
            --ts_filter_level 99.9 \
            -tranchesFile {tranches}  \
            -recalFile {recal} \
            -mode INDEL \
            -nt {num_jobs} \
            -o {output}".format(
                java=java,
                gatk=gatk,
                reference=reference,
                vcf=vcf,
                tranches=tranches,
                recal=recal,
                output=output,
                num_jobs=options.jobs
            ))

def mark_filtered_variants(vcf, output):
    """filter vcf"""
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
                output=output,
                input=vcf,
                reference=reference
            ))
# 
def remove_filtered(vcf,output):
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
                input=vcf,
                output=output
            ))

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
    files = glob.glob(options.bam_dir + '/*.bam')
    return [ os.path.splitext(os.path.basename(file))[0] for file in files ]

def generate_parameters():
    files = glob.glob(options.bam_dir + '/*.bam')
    parameters = []
    for file in files:
        prefix = os.path.splitext(os.path.basename(file))[0]
        parameters.append([None,prefix + '/' + prefix + '.bam', [prefix,os.path.abspath(file)]])
    for job_parameters in parameters:
            yield job_parameters

def get_num_files():
    files = glob.glob(options.bam_dir + '/*.bam')
    return len(files)

@files(generate_parameters)
def link(none, bam, extra):
    """Make working directory and make symlink to bam file"""
    if not os.path.exists(extra[0]):
        os.mkdir(extra[0])
    if not os.path.exists(bam):
        os.symlink(extra[1], bam) 

#@follows(link)
@transform(link, suffix(".bam"), '.bam.bai')
def index(input, output):
    """create bam index"""
    index_bam(input)

#@follows(index)
@transform(link, suffix(".bam"), '.dedup.bam')
def remove_dups(input, output):
    """Mark dups"""
    dup_removal_samtools(input, output)
    # remove(input)

#@follows(remove_dups)
@transform(remove_dups, suffix(".dedup.bam"), '.dedup.bam.bai')
def index_dups(input, output):
    """create bam index"""
    index_bam(input)

#@follows(index_dups)
@transform(index_dups, suffix(".dedup.bam.bai"), '.realign.intervals', r'\1.dedup.bam')
def find_realignment_intervals(foo, output, bam):
   """Find regions to be re-aligned due to indels"""
   find_realigns(bam, output)

#@follows(find_realignment_intervals)
@transform(find_realignment_intervals, suffix(".realign.intervals"), '.realigned.bam', r'\1.dedup.bam')
def indel_realigner(input, output, bam):
   """Re-aligns regions around indels"""
   run_realigner(bam,input,output)
   # remove(input)

#@follows(indel_realigner)
@transform(indel_realigner, suffix('.realigned.bam'), '.gatk.bam.recal_data.grp')
def recalibrate_baseq1(input, output):
    """Base quality score recalibration in bam file
        Part 1: count covariates"""
    index_bam(input)
    base_recalibrator(input, output)

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
def recalibrate_baseq2(inputs, output):
    """Base quality score recalibration in bam file
        Part 2: rewrite quality scores into a new bam file"""   
    print_recalibrated(inputs[0], inputs[1], output)
    # remove(inputs[0])


#
#
# QC measurements
#

#@follows(recalibrate_baseq2)
@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.quality_score')
def metric_quality_score_distribution(input,output):
    """docstring for metrics1"""
    bam_quality_score_distribution(input, output, output + '.pdf')

#@follows(recalibrate_baseq2)
@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.metrics')
def metric_alignment(input,output):
    """docstring for metrics1"""
    bam_alignment_metrics(input, output)

#@follows(recalibrate_baseq2)
@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.coverage.sample_summary', r'\1.coverage')
def metric_coverage(input, output, output_format):
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
                output=output_format,
                input=input,
                capture=options.capture_regions
            ))

# long running, memory demanding and not very useful
@merge(recalibrate_baseq2, 'multisample.coverage')
def metric_coverage_multisample(bams, output):
    bam_coverage_multisample_statistics(bams, output)
    cmd = "{java} -Xmx32g -jar {gatk} \
            -R {reference} \
            -T DepthOfCoverage \
            -o {output} \
            -L {capture} \
            -ct 8 -ct 20 -ct 30 \
            --omitDepthOutputAtEachBase --omitLocusTable \
        ".format(java=java, gatk=gatk,
                reference=reference,
                output=output,
                capture=options.capture_regions)

    for bam in bams:
        cmd += " -I {}".format(bam)

    run_cmd(cmd)


#
#
# bam reduce & variant calling
#


@jobs_limit(6)
#@follows(recalibrate_baseq2)
@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.reduced.bam')
def reduce_bam(input, output):
    reduce_reads(input,output)

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
# def split_snps_and_indels(input,output):
#     """Separates the vcf file into only snp calls"""
#     split_snps(input,output[0])
#     split_indels(input,output[1])

@follows(call_variants)
@files('multisample.gatk.vcf', ['multisample.gatk.snp.model','multisample.gatk.snp.model.tranches'])
def find_snp_tranches_for_recalibration(input,output):
    """Run variantRecalibration for trusted sites"""
    recalibrate_snps(input,output[0])

@follows(find_snp_tranches_for_recalibration)
@files('multisample.gatk.vcf', ['multisample.gatk.indel.model','multisample.gatk.indel.model.tranches'])
def find_indel_tranches_for_recalibration(input,output):
    """Run variantRecalibration for trusted sites"""
    recalibrate_indels(input,output[0])

@follows(find_indel_tranches_for_recalibration)
@files(['multisample.gatk.vcf','multisample.gatk.snp.model','multisample.gatk.snp.model.tranches'],'multisample.gatk.recalibratedSNPS.rawIndels.vcf')
def apply_recalibration_filter_snps(input,output):
    apply_recalibration_snps(input[0],input[1],input[2],output)
    # remove(input[0])
    #     remove(input[1])
    #     remove(input[2])

@follows(apply_recalibration_filter_snps)
@files(['multisample.gatk.recalibratedSNPS.rawIndels.vcf','multisample.gatk.indel.model','multisample.gatk.indel.model.tranches'],'multisample.gatk.preHardFiltering.vcf')
def apply_recalibration_filter_indels(input,output):
    apply_recalibration_indels(input[0],input[1],input[2],output)

@follows(apply_recalibration_filter_indels)
@files('multisample.gatk.preHardFiltering.vcf', 'multisample.gatk.markedHardFiltering.vcf')
def apply_filter(input,output):
    """docstring for apply_indel_filter"""
    mark_filtered_variants(input,output)
    # remove(input)

@follows(apply_filter)
@files('multisample.gatk.markedHardFiltering.vcf', 'multisample.gatk.analysisReady.vcf')
def filter_variants(infile,outfile):
    """Filter variants that passed"""
    remove_filtered(infile,outfile)
    # remove(infile)

@follows(filter_variants)
@transform(filter_variants,suffix('.analysisReady.vcf'),'.analysisReady.exome.vcf')
def final_calls(input,output):
    """Produce the final variant calls"""
    output = output[:-10]
    filter_by_exome_region(input, output)
    rename('multisample.gatk.analysisReady.recode.vcf','multisample.gatk.analysisReady.exome.vcf')

def split_snp_parameters():
    exome_vcf = 'multisample.gatk.analysisReady.exome.vcf'
    for id in get_sample_ids():
        yield [exome_vcf, id + '/' + id + '.exome.vcf', id]

@posttask(cleanup_files)
@follows('final_calls')
@files(split_snp_parameters)
def split_snps(input, output, sample):
    """ Split variants by sample, and use sample-specific statistics to filter: AD and DP"""
    split_variants_by_sample(input, sample, output, ad_threshold=5, dp_threshold=8)



###########
# Annovar annotation / filtering
###########


# reference dbs
annovar_human_db = os.path.join(script_path,'../tools/annovar/humandb')
annovar_1000genomes_eur = '1000g2012apr_eur'
annovar_inhouse_db = 'common_inhouse_variants_jan2014.txt'

@transform(split_snps, formatter('.*/(?P<SAMPLE_ID>[^/]+).exome.vcf'),
			'{subpath[0][1]}/annotated-with-annovar/{SAMPLE_ID[0]}.avinput')
def prepare_annovar_inputs(vcf, output):
    try: os.mkdir('annotated-with-annovar')
    except (OSError): pass # dir exists
 
    run_cmd("convert2annovar.pl {vcf} -format vcf4 -withzyg -includeinfo -outfile {out} \
        ".format(vcf=vcf, out=output))
    

@split(final_calls,'annotated-with-annovar/*.avinput') 
def prepare_annovar_inputs2(multisample_vcf, outputs):
    """ create an annovar file for every sample. needs to be run separately from the rest"""
    os.mkdir('annotated-with-annovar')
    output_prefix = 'annotated-with-annovar/sample'
    run_cmd("convert2annovar.pl {vcf} -format vcf4 -withzyg -includeinfo \
        -genoqual 3 -coverage 6 -allsample -outfile {output_prefix}".format(
        vcf=multisample_vcf, 
        output_prefix=output_prefix))
    
 
def annotate_variants_with_functional_change(input_file, output_prefix):
    run_cmd("annotate_variation.pl -buildver hg19 -outfile {outfile_prefix} {input_file} {annodb}".format(
        outfile_prefix=output_prefix, 
        input_file=input_file, 
        annodb=annovar_human_db))


@transform(prepare_annovar_inputs, suffix('.avinput'), ['.avinput.variant_function.stats','.avinput.exonic_variant_function.stats'])
def annotate_function_of_raw_variants(input, outputs):
    """ annotate functional change in raw variants """
    annotate_variants_with_functional_change(input_file=input, output_prefix=input)
    # calculate stats on files created by annovar - output files without ".stats" suffix
    run_cmd("cut -f 1 {f} | sort | uniq -c > {f}.stats".format(f=outputs[0][:-len('.stats')]))
    run_cmd("cut -f 2 {f} | sort | uniq -c > {f}.stats".format(f=outputs[1][:-len('.stats')]))
    # remove the annovar files
    remove(outputs[0][:-len('.stats')])
    remove(outputs[1][:-len('.stats')])


@transform(prepare_annovar_inputs, suffix('.avinput'), 
					['.avinput.common_inhouse_filtered', 
					 '.avinput.common_inhouse_dropped'])
def filter_common_inhouse(input, outputs):
    """ filter variants found in the inhouse database. OBS! output specifies the filename after rename """    
    run_cmd("annotate_variation.pl -build hg19 -filter -dbtype generic -genericdbfile {inhouse} \
        -outfile {outfile} {input} {annodb}".format(
        inhouse=annovar_inhouse_db, 
        outfile=input, 
        input=input, 
        annodb=annovar_human_db))

    for output_file in outputs: 
	       rename(output_file.replace('common_inhouse','hg19_generic'), output_file)

    
@transform(filter_common_inhouse, suffix('.common_inhouse_filtered'),
					['.common_inhouse_filtered.hg19_EUR.sites.2012_04_filtered',
					 '.common_inhouse_filtered.hg19_EUR.sites.2012_04_dropped'])
def filter_common_1000genomes(inputs, outputs):
    """ filter common 1000 genomes variants """
    filtered = inputs[0]                      # take only the filtered file, leave dropped
    run_cmd("annotate_variation.pl -build hg19 -filter -dbtype {eur1kg} \
        -maf {maf} -outfile {output_prefix} {input_file} {annodb}".format(
        eur1kg=annovar_1000genomes_eur,
        maf=0.005, 
        output_prefix=filtered,
        input_file=filtered, 
        annodb=annovar_human_db))


@transform(filter_common_1000genomes, suffix('.hg19_EUR.sites.2012_04_filtered'), 
	   				['.hg19_EUR.sites.2012_04_filtered.variant_function',
					 '.hg19_EUR.sites.2012_04_filtered.exonic_variant_function',
					 '.hg19_EUR.sites.2012_04_filtered.variant_function.stats',
					 '.hg19_EUR.sites.2012_04_filtered.exonic_variant_function.stats'])
def annotate_function_of_rare_variants(inputs, outputs):
    """ annotate functional change in rare variants """
    filtered = inputs[0] 			# use only the filtered input file, leave dropped
    annotate_variants_with_functional_change(input_file=filtered, output_prefix=filtered)
    # calculate stats on files created by annovar
    run_cmd("cut -f 1 {f} | sort | uniq -c > {f}.stats".format(f=outputs[0]))
    run_cmd("cut -f 2 {f} | sort | uniq -c > {f}.stats".format(f=outputs[1]))


@transform(annotate_function_of_rare_variants, 
           formatter(".*/(?P<SAMPLE_ID>[^/]+).avinput.common_inhouse_filtered.hg19_EUR.sites.2012_04_filtered.variant_function", None, None, None),
           ['{path[0]}/annotated-tables/{SAMPLE_ID[0]}.rare_coding_and_splicing.avinput', 
            '{path[0]}/annotated-tables/{SAMPLE_ID[0]}.rare_coding_and_splicing.avinput.hg19_multianno.csv'])
def produce_variant_annotation_table(inputs, outputs):
    """ produce a table of various annotations per variant """
	
    dir = 'annotated-with-annovar/annotated-tables'
    try:
	       os.mkdir(dir)
    except (OSError):
	       pass # dir exists
    
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
        if lsplit[0] == 'splicing':
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

omim_gene_phenotype_map = get_omim_gene_phenotype_map(os.path.join(script_path,'../tools/annovar/omim/genemap2.txt'))


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
		   '20','21','22','X','Y','GL000209.1']

    hetz = open(table_files[0], 'w')
    homz = open(table_files[1], 'w')
    hetz.write('sample\t'+string.join(chromosomes,'\t')+'\n')
    homz.write('sample\t'+string.join(chromosomes,'\t')+'\n')
    for fname in infiles:
	sample_id=os.path.basename(fname).split('.')[1]
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

@merge([annotate_function_of_raw_variants, annotate_function_of_rare_variants], 'all_samples_exonic_variant_stats.tsv')
def produce_variant_stats_table(infiles, table_file):
    """ produce a table of per-sample counts of different type of exonic variants """
    # split the input files per task
    sample_no = len(infiles)/2
    raw_variant_files = infiles[0:sample_no]
    rare_variant_files = infiles[sample_no:]

    out = open(table_file,'w')    
    out.write('sample\traw_exonic\trare_exonic\traw_synonymous\trare_synonymous\n')
    for i in range(0,sample_no):
        out.write(os.path.basename(raw_variant_files[i][0]).split('.')[1])
        for fname in [raw_variant_files[i][0], rare_variant_files[i][2]]: # exonic variant stats of raw variants and rare variants
            f=open(fname)	
            for l in f.xreadlines():
                if l.find("exonic")>0:
                    out.write('\t'+l.split()[0])
                    break
            f.close()
        for fname in [raw_variant_files[i][1], rare_variant_files[i][3]]: # synonymous variants stats of raw and rare variants
            f=open(fname)
            for l in f.xreadlines():
                if l.find(" synonymous")>0:
                    out.write('\t'+l.split()[0])
                    break            
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

