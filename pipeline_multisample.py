#!/usr/bin/env python
"""

    pipeline_multisample.py
                        [--settings PATH]
                        [--log_file PATH]
                        [--verbose]
                        [--target_tasks]
                        [--jobs N]
                        [--just_print]
                        [--flowchart]
                        [--key_legend_in_graph]
                        [--forced_tasks]
                        [--run_on_bcl_tile TILE_REGEX]

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

    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --settings PIPELINE_SETTINGS.CFG --target_task TASK [more_options]")
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
    parser.add_option("--run_on_bcl_tile", dest="run_on_bcl_tile",
                        type="string", 
                        default=None,
                        help="Use only this tile when doing bcl2fastq conversion. For testing purposes.")
    

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

    # Root dirs
    data_root = config.get('Docker','data-root')
    reference_root = config.get('Docker','reference-root')
      
    # optional results and fastq archive dirs  
    results_archive = None
    try:
        results_archive = config.get('Docker','results-archive')
    except ConfigParser.NoOptionError:
        print 'No results-archive provided. Results will not be archived outside of the current (working) directory.'
    
    fastq_archive = None
    try:
        fastq_archive = config.get('Docker','fastq-archive')
    except ConfigParser.NoOptionError:
        print 'No fastq-archive provided. Fastq files will not be archived outside of the current (working) directory.'

    
    
    tmp_dir = None
    try:
        tmp_dir = config.get('Docker','tmp-dir')
    except ConfigParser.NoOptionError:
        print 'No tmp-dir provided. Container\'s /tmp will be used.'
    
    
    # Docker executable and args
    docker_bin = config.get('Docker','docker-binary')
    docker_args = config.get('Docker', 'docker-args')
    docker_args += " -v " + ":".join([data_root,data_root,"ro"])
    docker_args += " -v " + ":".join([reference_root,reference_root,"ro"])
    if results_archive != None:
        docker_args += " -v " + ":".join([results_archive,results_archive,"rw"])
        
    if tmp_dir != None: 
        docker_args += " -v " + ":".join([tmp_dir,tmp_dir,"rw"])
    else: # set the default value if the tmp-dir was unset
        tmp_dir = "/tmp"
        
    cwd = os.getcwd()
    docker_args += " -v " + ":".join([cwd,cwd,"rw"])
    docker_args += " -w " + cwd
    docker = " ".join([docker_bin, docker_args]) 
    
    # Inputs 
    input_bams = os.path.join(data_root, config.get('Inputs','input-bams'))
    input_dir  = os.path.join(data_root, config.get('Inputs','input-dir'))
    
    # variant calls from other projects to call together with (semicolon separated list)
    try: 
        call_with_gvcfs = [ os.path.join(reference_root, path.strip()) for path in config.get('Inputs','call-with-gvcfs').split(";") ]
        #print 'Calling will be done together with:'
        #for p in call_with_gvcfs:
        #	print '\t',p
    except (ConfigParser.NoOptionError): 
        call_with_gvcfs = [] 
    
    # reference files
    reference = os.path.join(reference_root, config.get('Resources','reference-genome'))
    dbsnp = os.path.join(reference_root, config.get('Resources','dbsnp-vcf'))
    hapmap = os.path.join(reference_root, config.get('Resources','hapmap-vcf'))
    omni = os.path.join(reference_root, config.get('Resources','1000genomes-omni-vcf'))
    snps_1kg = os.path.join(reference_root, config.get('Resources','1000genomes-snps-vcf'))
    indels_1kg = os.path.join(reference_root, config.get('Resources','1000genomes-indels-vcf'))
    mills = os.path.join(reference_root, config.get('Resources','mills-indels-vcf'))
    capture = os.path.join(reference_root, config.get('Resources','capture-regions-bed'))
    capture_qualimap = os.path.join(reference_root, config.get('Resources','capture-regions-bed-for-qualimap'))
    exome = os.path.join(reference_root, config.get('Resources', 'exome-regions-bed'))
    gene_coordinates = os.path.join(reference_root, config.get('Resources', 'gene-coordinates'))
    
    adapters = os.path.join(reference_root, config.get('Resources', 'adapters-fasta'))
    
    # tools
    bcl2fastq = config.get('Tools','bcl2fastq')
    trimmomatic = config.get('Tools', 'trimmomatic') 
    bwa = config.get('Tools','bwa')
    samtools = config.get('Tools','samtools')
    picard = config.get('Tools','picard-tools')
    qualimap = config.get('Tools','qualimap')
    gatk = config.get('Tools','gatk')
    vcftools = config.get('Tools','vcftools')


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#  Common functions 


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

import drmaa
drmaa_session = drmaa.Session()
drmaa_session.initialize()

from ruffus.drmaa_wrapper import run_job, error_drmaa_job

def _run_cmd(cmd, cpus, mem_per_cpu, run_locally):
    stdout, stderr = '', ''
    job_options = '--ntasks=1 \
                    --cpus-per-task={cpus} \
                    --mem-per-cpu={mem} \
                    '.format(cpus=cpus, mem=mem_per_cpu)
                   
    try:
        stdout, stderr = run_job(cmd, 
                                 job_other_options=job_options,
                                 run_locally = run_locally, 
                                 retain_job_scripts = True, job_script_directory = 'drmaa/',
                                 drmaa_session = drmaa_session)
    except error_drmaa_job as err:
        raise Exception("\n".join(map(str, ["Failed to run:", cmd, err, stdout, stderr])))
    
def run_cmd(cmd, args, interpreter_args=None, cpus=1, mem_per_cpu=1024, run_locally=False, dockerize=True):
    
    if dockerize:

#        full_cmd = "{docker} {cmd} \"{args}\"".format(docker=docker,cmd=cmd,args=args)
#        if interpreter_args!=None and interpreter_args.strip()!="":
#            full_cmd += " \"{}\"".format(interpreter_args)
            
        if interpreter_args==None:
            full_cmd = "{docker} {cmd} {args}".format(docker=docker,cmd=cmd,args=args)
        else:
            full_cmd = "{docker} {cmd} \"{args}\" \"{iargs}\" \
                    ".format(docker=docker, cmd=cmd, args=args, iargs=interpreter_args)
    else: 
        
        if interpreter_args==None or interpreter_args.strip()=="":
            full_cmd = "{cmd} {args}".format(cmd=cmd, args=args)
        else:
            raise Exception('Not implemented interpreter args in not-dockerized execution.\nWhy would you want to do it anyway?')
            #full_cmd = "{cmd} \"{args}\" \"{iargs}\"".format(cmd=cmd, args=args, iargs=interpreter_args)

    _run_cmd(full_cmd, cpus, mem_per_cpu, run_locally)

def rename(old_file, new_file):
    """rename file"""
    os.rename(old_file, new_file)

def remove(f):
    """remove file"""
    os.remove(f)
                            
def index_bam(bam):
    """Use samtools to create an index for the bam file"""
    run_cmd(samtools, "index %s" % bam)
                          
def bam_quality_score_distribution(bam,qs,pdf):
    """Calculates quality score distribution histograms"""
    run_cmd(picard, "QualityScoreDistribution \
                    CHART_OUTPUT={chart} \
                    OUTPUT={out} \
                    INPUT={bam} \
                    VALIDATION_STRINGENCY=SILENT \
                    ".format(chart=pdf, out=qs, bam=bam),
            interpreter_args="")

def bam_alignment_metrics(bam,metrics):
    """Collects alignment metrics for a bam file"""
    run_cmd(picard, "CollectAlignmentSummaryMetrics \
                    REFERENCE_SEQUENCE={ref} \
                    OUTPUT={out} \
                    INPUT={bam} \
                    VALIDATION_STRINGENCY=SILENT \
                    ".format(ref=reference, out=metrics, bam=bam),
            interpreter_args="")
    
def bam_target_coverage_metrics(input_bam, output):
    """ Calculates and outputs bam coverage statistics """
    run_cmd(gatk, "-R {reference} \
                    -T DepthOfCoverage \
                    -o {output} \
                    -I {input} \
                    -L {capture} \
                    -ct 8 -ct 20 -ct 30 \
                    --omitDepthOutputAtEachBase --omitLocusTable \
                    ".format(reference=reference,
                             output=output,
                             input=input_bam,
                             capture=capture), 
            interpreter_args="-Xmx4g")

def bam_gene_coverage_metrics(input_bam, output):
    """ Calculates and outputs bam coverage statistics """
    run_cmd(gatk, "-R {reference} \
                    -T DepthOfCoverage \
                    -o {output} \
                    -I {input} \
                    -L {capture} \
                    -geneList {genes} \
                    -ct 5 -ct 10 -ct 20 \
                    --omitDepthOutputAtEachBase --omitLocusTable \
                    ".format(reference=reference,
                             output=output,
                             input=input_bam,
                             capture=capture,
                             genes=gene_coordinates), 
            interpreter_args="-Xmx4g")


def qualimap_bam(input_bam, output_dir):
    """ Generates Qualimap bam QC report """
    # create necessary folders first
    if not os.path.exists('qc'): os.mkdir('qc')
    if not os.path.exists('qc/qualimap/'): os.mkdir('qc/qualimap')
    if not os.path.exists(output_dir): os.mkdir(output_dir)
    run_cmd(qualimap, "bamqc -bam {bam} \
                        -c -outformat PDF \
                        -gff {target} \
                        -gd HUMAN -os \
                        -outdir {dir} &> {dir}/qualimap.err \
                        ".format(bam=input_bam,
                                target=capture_qualimap,
                                dir=output_dir),
            interpreter_args="")

    

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

from ruffus import *



def get_sample_ids():
    #files = glob.glob(input_bams)
    #return [ os.path.splitext(os.path.basename(f))[0] for f in files ]
    
    files = glob.glob(os.path.join(cwd,'*','*.gvcf'))
    return [ os.path.splitext(os.path.basename(f))[0] for f in files ]
    # expect files in form /path/to/data/dir/[SAMPLE_ID]_[LANE_ID]_R[1|2]_001.bcl
    # and return a list of unique sorted SAMPLE_IDs
   # ids = [ (os.path.splitext(os.path.basename(f))[0]).split('_')[0] for f in files ]
   # uniq_ids = list(set(ids))
   # uniq_ids.sort()
   # return uniq_ids

def get_num_files():
    return len(get_sample_ids())

def are_fastqs_converted(_,__):
#    fastqs = glob.glob(os.path.join(cwd,'fastqs','*.fastq.gz'))
#    print(fastqs)

#    for fastq in fastqs:
#        if not fastq.startswith('Undetermined'): 
#            return False, 'FASTQ file %s found. Skipping bcl2fastq conversion' % os.path.basename(fastq)
#   return True, 'No FASTQ files found'
    if os.path.exists(os.path.join(cwd,'fastqs','completed')):
        return False, 'Found bcl2fastq completion flag'
    return True, 'Missing bcl2fastq completion flag'


#
# Preprocess the reads, only if fastqs don't exist
#
@files(input_dir, os.path.join(cwd,'fastqs','completed'))
#@check_if_uptodate(are_fastqs_converted)
@posttask(touch_file(os.path.join(cwd,'fastqs','completed')))
def bcl2fastq_conversion(run_directory):
    """ Run bcl2fastq conversion and create fastq files in the run directory"""
    out_dir = os.path.join(cwd,'fastqs')
    interop_dir = os.path.join(out_dir,'InterOp')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir) 

    # r, w, d, and p specify numbers of threads to be used for each of the concurrent subtasks of the conversion (see bcl2fastq manual) 
    args = "-R {indir} -o {outdir} --interop-dir={interopdir} -r1 -w1 -d2 -p4 \
            ".format(indir=run_directory, outdir=out_dir, interopdir=interop_dir)
    if options.run_on_bcl_tile != None:
        args += " --tiles %s" % options.run_on_bcl_tile
    run_cmd(bcl2fastq, args, cpus=8, mem_per_cpu=2048)
    
    # touch a flag indicating that it has finished conversion
    #open(os.path.join(out_dir,'completed'),'w').close()

@active_if(fastq_archive != None)
#@transform(bcl2fastq_conversion, formatter(".+/?P<RUN_ID>/fastqs/completed"), os.path.join(fastq_archive,'{RUN_ID[0]}'))
@transform(bcl2fastq_conversion, formatter(".+/(?P<RUN_ID>[^/]+)/fastqs/completed"), str(fastq_archive) + '/{RUN_ID[0]}')
def archive_fastqs(input, archive_dir):
    """ Archive fastqs """
    # if optional fastq-archive-root was not provided - do nothing
    #if fastq_archive == None: return
    
    fq_dir = os.path.dirname(input)

#    run_name = os.path.basename(cwd)
#    arch_path = os.path.join(fastq_archive, run_name)
#    if not os.path.exists(arch_path):
#        os.mkdir(arch_path)

    import shutil
    shutil.move(fq_dir, archive_dir)
    os.mkdir(fq_dir)
    for f in glob.glob(os.path.join(archive_dir,'*.fastq.gz')):
        os.symlink(f, os.path.join(fq_dir,os.path.basename(f)))


#
# Prepare directory structure move in the input fastq files
# Expected format:
#    /path/to/file/[SAMPLE_ID]_[anything_except_path_delimiter].fastq.gz
#
@jobs_limit(1)    # to avoid problems with simultanous creation of the same sample dir
@follows(archive_fastqs)
@subdivide(os.path.join(cwd,'fastqs','*.fastq.gz'),
           formatter('.+/(?P<SAMPLE_ID>[^_/]+)_S[1-9][0-9]?_L\d\d\d_R[12]_001\.fastq\.gz$'), 
           '{SAMPLE_ID[0]}/{basename[0]}{ext[0]}')
def link_fastqs(fastq_in, fastq_out):
    """Make working directory for every sample and move fastq files in"""
    if not os.path.exists(os.path.dirname(fastq_out)):
        os.mkdir(os.path.dirname(fastq_out))
    if not os.path.exists(fastq_out):
        os.symlink(fastq_in, fastq_out) 

    
    
#
# Input FASTQ filenames are expected to have following format:
#    [SAMPLE_ID]_[S_NUM]_[LANE_ID]_[R1|R2]_001.fastq.gz
# In this step, the two FASTQ files matching on the [SAMPLE_ID]_[S_ID]_[LANE_ID] will be trimmed together (R1 and R2). 
# The output will be written to two FASTQ files
#    [SAMPLE_ID]_[LANE_ID].fq1.gz
#    [SAMPLE_ID]_[LANE_ID].fq2.gz
#
@collate(link_fastqs, regex(r'([^_]+)_S[1-9]\d?_(L\d\d\d)_R[12]_001\.fastq\.gz$'),  r'\1_\2.fq1.gz')
def trim_reads(inputs, output):
    outfq1 = output
    outfq2 = output.replace('fq1.gz','fq2.gz')
    unpaired = [outfq1.replace('fq1.gz','fq1_unpaired.gz'), outfq2.replace('fq2.gz','fq2_unpaired.gz')]               
    # logfile = output.replace('fq1.gz','trimmomatic.log')
    # -trimlog {log} \
    # log=logfile
    args = "PE -phred33 -threads 1 \
            {in1} {in2} {out1} {unpaired1} {out2} {unpaired2} \
            MINLEN:36 \
            ILLUMINACLIP:{adapter}:2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15".format(in1=inputs[0], in2=inputs[1],
                                       out1=outfq1, out2=outfq2,
                                       unpaired1=unpaired[0], unpaired2=unpaired[1],
                                       adapter=adapters)
    max_mem = 2048
    run_cmd(trimmomatic, args, interpreter_args="-Xmx"+str(max_mem)+"m", cpus=1, mem_per_cpu=max_mem)


#
# FASTQ filenames are expected to have following format:
#    [SAMPLE_ID]_[LANE_ID].fq[1|2].gz
# In this step, the fq1 file coming from trim_reads is matched with the fq2 file and mapped together. 
# The output will be written to SAM file:
#    [SAMPLE_ID]_[LANE_ID].sam
#
#@collate(trim_reads, regex(r"([^_]+_[^_]+)\.fq[12]\.gz$"), r'\1.sam')
@transform(trim_reads, suffix('.fq1.gz'), add_inputs(r'\1.fq2.gz'), '.sam')
def align_reads(fastqs, sam):
    threads = 2
    
    # construct read group information from fastq file name (assuming [SAMPLE_ID]_[LANE_ID].fq[1|2].gz format)
    sample_lane = os.path.basename(fastqs[0])[0:-len(".fq1.gz")]
    sample = "_".join(sample_lane.split("_")[0:-1])
    lane = sample_lane.split("_")[-1]
    read_group = "@RG\\tID:{id}\\tSM:{sm}\\tLB:{lb}\\tPL:{pl}\\tPU:{pu} \
                 ".format(id=sample_lane, sm=sample, lb=sample, pl="ILLUMINA", pu=lane)
                 
    args = "mem -t {threads} -R '{rg}' {ref} {fq1} {fq2} > {sam} \
            ".format(threads=threads, rg=read_group, ref=reference, 
                     fq1=fastqs[0], fq2=fastqs[1], sam=sam)
    run_cmd(bwa, args, cpus=threads, mem_per_cpu=8192/threads)
    

#
# SAM filenames are expected to have following format:
#    [SAMPLE_ID]_[LANE_ID].sam
# In this step, all SAM files matching on the SAMPLE_ID will be merged into one BAM file:
#    [SAMPLE_ID].bam
#
@collate(align_reads, regex(r"([^_]+).+\.sam$"),  r'\1.bam')
def merge_lanes(lane_sams, bam):
    args = "MergeSamFiles O={bam} \
            ASSUME_SORTED=false \
            ".format(bam=bam)
    # include all sam files as args
    for sam in lane_sams:
        args += " I={sam}".format(sam=sam)
        
    run_cmd(picard, args, interpreter_args="-Xmx8g", cpus=4, mem_per_cpu=2048)
    

#
# For bam-level entry to the pipeline. 
#
def generate_bam_inputs():
    files = glob.glob(input_bams)
    parameters = []
    for f in files:
        prefix = os.path.splitext(os.path.basename(f))[0]
        parameters.append([None,prefix + '/' + prefix + '.bam', [prefix,os.path.abspath(f)]])
    for job_parameters in parameters:
            yield job_parameters


def clean_fastqs_and_sam():
    """ Remove the trimmed fastq files, and SAM files. Links to original fastqs are kept """
    for f in glob.glob(os.path.join(cwd,'*','*.fq[12]*.gz')):
        os.remove(f)
    for f in glob.glob(os.path.join(cwd,'*','*.sam')):
        os.remove(f)


@posttask(clean_fastqs_and_sam)
@transform(merge_lanes, suffix(".bam"), '.bam.bai')
def index(bam, output):
    """Create raw bam index"""
    index_bam(bam)



#
#
# QC the raw bam files
#

@follows(index)
@transform(merge_lanes, suffix(".bam"), '.quality_score')
def qc_raw_bam_quality_score_distribution(input_bam, output):
    """docstring for metrics1"""
    bam_quality_score_distribution(input_bam, output, output + '.pdf')

@follows(index)
@transform(merge_lanes, suffix(".bam"), '.metrics')
def qc_raw_bam_alignment_metrics(input_bam, output):
    """docstring for metrics1"""
    bam_alignment_metrics(input_bam, output)
    
@follows(index)
@transform(merge_lanes, suffix(".bam"), '.target_coverage.sample_summary', r'\1.target_coverage')
def qc_raw_bam_target_coverage_metrics(input_bam, output, output_format):
    bam_target_coverage_metrics(input_bam, output_format)

@follows(index)
@transform(merge_lanes, suffix('.bam'), '.gene_coverage.sample_summary', r'\1.gene_coverage')
def qc_raw_bam_gene_coverage_metrics(input_bam, output, output_format):
    bam_gene_coverage_metrics(input_bam, output_format)

@follows(index)
@transform(merge_lanes, formatter(".*/(?P<SAMPLE_ID>[^/]+).bam"), '{subpath[0][1]}/qc/qualimap/{SAMPLE_ID[0]}')
def qc_raw_bam_qualimap_report(input_bam, output_dir):
    qualimap_bam(input_bam, output_dir)

@follows(qc_raw_bam_target_coverage_metrics, qc_raw_bam_qualimap_report)
def raw_bam_qc():
    """ Aggregates raw bam quality control steps """
    pass





#
#
# Prepare the bam files for variant calling
#

def dup_removal_picard(bam,output):
    """Use Picard to remove duplicates"""
    args = "CleanSam \
           INPUT={bam} \
           OUTPUT={out} \
           VALIDATION_STRINGENCY=LENIENT \
           VERBOSITY=ERROR CREATE_INDEX=TRUE \
           ".format(bam=bam, out=output)
    run_cmd(picard, args, interpreter_args="-Xmx4g", mem_per_cpu=4096)
   
def dup_mark_picard(bam,output):
    """Use Picard to mark duplicates"""
    args = "MarkDuplicates \
            TMP_DIR={tmp} \
            REMOVE_DUPLICATES=true \
            INPUT={bam} \
            OUTPUT={out} \
            METRICS_FILE={bam}.dup_metrics \
            VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR CREATE_INDEX=true \
            ".format(tmp=tmp_dir, 
                     bam=bam, 
                     out=output)
    run_cmd(picard, args, interpreter_args="-Xmx4g", mem_per_cpu=4096)


@follows(index)
@transform(merge_lanes, suffix(".bam"), '.dedup.bam')
def remove_dups(bam, output):
    """Mark duplicates"""
    #run_cmd(samtools,"rmdup %s %s" % (bam,output))
    dup_mark_picard(bam, output)
    # remove(input)


#@follows(remove_dups)
@transform(remove_dups, suffix(".dedup.bam"), '.dedup.bam.bai')
def index_dups(bam, output):
    """Create bam index"""
    index_bam(bam)


#@transform(index_dups, suffix(".dedup.bam.bai"), '.realign.intervals', r'\1.dedup.bam')
#def find_realignment_intervals(foo, intervals, input_bam):
@transform(remove_dups, suffix(".dedup.bam"), '.realign.intervals')
def find_realignment_intervals(input_bam, intervals):
    """ Generate realignment intervals """
    args = "-T RealignerTargetCreator \
            -I {bam} \
            -R {reference} \
            -known {indels1} \
            -known {indels2} \
            -o {out}".format(bam=input_bam, 
                             reference=reference, 
                             indels1=indels_1kg, 
                             indels2=mills, 
                             out=intervals)
    run_cmd(gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx4g" % tmp_dir, mem_per_cpu=4096)


@transform(find_realignment_intervals, suffix(".realign.intervals"), add_inputs(r'\1.dedup.bam'), '.realigned.bam')
def indel_realigner(inputs, realigned_bam):
    intervals_file = inputs[0]
    input_bam = inputs[1]
    """Re-aligns regions around indels"""
    args = "-T IndelRealigner \
            -I {bam} \
            -R {reference} \
            -targetIntervals {intervals} \
            -known {indels1} \
            -known {indels2} \
            -o {out}".format(bam=input_bam, 
                             reference=reference, 
                             intervals=intervals_file, 
                             indels1=indels_1kg, 
                             indels2=mills, 
                             out=realigned_bam)
    run_cmd(gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx4g" % tmp_dir, mem_per_cpu=4096)


@transform(indel_realigner, suffix('.realigned.bam'), '.gatk.bam.recal_data.grp')
def recalibrate_baseq1(input_bam, output):
    """Base quality score recalibration in bam file - first pass """
    #index_bam(input_bam)
    args = "-T BaseRecalibrator \
           -R {reference} \
           -knownSites {dbsnp} \
           -I {bam} \
           -o {out}".format(reference=reference, 
                            dbsnp=dbsnp, 
                            bam=input_bam, 
                            out=output)
    run_cmd(gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx6g" % tmp_dir, mem_per_cpu=6144)
    

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
    args = "-T PrintReads \
            -R {reference} \
            -I {bam} \
            --out {out} \
            -BQSR {recal}".format(reference=reference, 
                                  bam=bam, 
                                  out=output_bam, 
                                  recal=recal_data) 
    run_cmd(gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx4g" % tmp_dir, mem_per_cpu=4096)
    
    # remove(inputs[0])


#
#
# gatk.bam-level QC measurements
#

@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.quality_score')
def qc_gatk_bam_quality_score_distribution(input_bam, output):
    """docstring for metrics1"""
    bam_quality_score_distribution(input_bam, output, output + '.pdf')

@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.metrics')
def qc_gatk_bam_alignment_metrics(input_bam, output):
    """docstring for metrics1"""
    bam_alignment_metrics(input_bam, output)

@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.target_coverage.sample_summary', r'\1.target_coverage')
def qc_gatk_bam_target_coverage_metrics(input_bam, output, output_format):
    bam_target_coverage_metrics(input_bam, output_format)

@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.gene_coverage.sample_summary', r'\1.gene_coverage')
def qc_gatk_bam_gene_coverage_metrics(input_bam, output, output_format):
    bam_gene_coverage_metrics(input_bam, output_format)

@transform(recalibrate_baseq2, formatter(".*/(?P<SAMPLE_ID>[^/]+).bam"), '{subpath[0][1]}/qc/qualimap/{SAMPLE_ID[0]}')
def qc_gatk_bam_qualimap_report(input_bam, output_dir):
    qualimap_bam(input_bam, output_dir)

@follows(qc_gatk_bam_gene_coverage_metrics)
def gatk_bam_qc():
    """ Aggregates gatk_bam quality control steps """
    pass


#
#
# Genotyping using HaplotypeCaller
#

@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.gvcf')
def call_haplotypes(bam, output_gvcf):
    """Perform variant calling using GATK HaplotypeCaller"""
    args = "-T HaplotypeCaller \
            -R {reference} \
            -I {bam} \
            -o {gvcf} \
            --emitRefConfidence GVCF \
            --variant_index_type LINEAR \
            --variant_index_parameter 128000 \
            -minPruning 4 \
            -L {target} \
            --dbsnp {dbsnp} \
            ".format(reference=reference, 
                     bam=bam, 
                     gvcf=output_gvcf, 
                     target=exome, 
                     dbsnp=dbsnp)
#             -log {gvcf}.log \

    run_cmd(gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx6g" % tmp_dir, mem_per_cpu=6144)


@merge(call_haplotypes, 'multisample.gatk.gvcf')
def merge_gvcfs(gvcfs, merged_gvcf):
    """Combine the per-sample GVCF files into one project-wide GVCF""" 
    args = "-T CombineGVCFs \
            -R {reference} \
            -o {out} \
            -log {out}.log".format(reference=reference, out=merged_gvcf)
       
    for gvcf in gvcfs:
        args = args + " --variant {}".format(gvcf)
    
    run_cmd(gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx4g" % tmp_dir, mem_per_cpu=4096)
    

@originate(call_with_gvcfs)
def call_with_gvcfs_as_arguments_task(gvcfs):
    pass

@merge([merge_gvcfs, call_with_gvcfs_as_arguments_task], 'multisample.gatk.vcf')
def genotype_gvcfs(gvcfs, output):
    """ Genotype this project's merged GVCF together with other project-wide GVCF files (provided in settings) """
    args = "-T GenotypeGVCFs \
            -R {reference} \
            -o {out} \
            -log {out}.log".format(reference=reference, out=output)

    # if there are any external gvcfs to call with, include them
    for gvcf in gvcfs:
        args = args + " --variant {}".format(gvcf)

    run_cmd(gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx8g" % tmp_dir, mem_per_cpu=8192)



#
#
# Variant recalisbration and filtering
# 

#@follows(call_variants)
@follows(genotype_gvcfs)
@files('multisample.gatk.vcf', ['multisample.gatk.snp.model','multisample.gatk.snp.model.tranches'])
def find_snp_tranches_for_recalibration(vcf,outputs):
    """Runs VariantRecalibrator on snps"""
    args = "-T VariantRecalibrator \
            -R {reference}  \
            -input {input} \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} \
            -resource:omni,known=false,training=true,truth=true,prior=12.0 {omni} \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 {snps_1kg} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an DP \
            -mode SNP \
            -recalFile {output} \
            -tranchesFile {tranches} \
            -rscriptFile {plots}\
            -nt {num_jobs}".format(
                reference=reference,
                hapmap=hapmap,
                omni=omni,
                snps_1kg=snps_1kg,
                dbsnp=dbsnp,
                input=vcf,
                output=outputs[0],
                tranches=outputs[1],
                plots=outputs[0]+'.plots.R',
                num_jobs=options.jobs)
            
    if get_num_files() > 10:
        args += " -an InbreedingCoeff"
        
    run_cmd(gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx16g" % tmp_dir, cpus=options.jobs, mem_per_cpu=16386/options.jobs)


@follows(find_snp_tranches_for_recalibration)
@files('multisample.gatk.vcf', ['multisample.gatk.indel.model','multisample.gatk.indel.model.tranches'])
def find_indel_tranches_for_recalibration(vcf,outputs):
    """Runs VariantRecalibrator on indels"""
    args = "-T VariantRecalibrator \
            -R {reference}  \
            -input {input} \
            --maxGaussians 4 \
            -resource:mills,known=false,training=true,truth=true,prior=12.0 {mills} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} \
            -an QD -an DP -an FS -an ReadPosRankSum -an MQRankSum \
            -mode INDEL \
            -recalFile {output} \
            -tranchesFile {tranches} \
            -rscriptFile {plots}\
            -nt {num_jobs}".format(
                reference=reference,
                mills=mills,
                dbsnp=dbsnp,
                input=vcf,
                output=outputs[0],
                tranches=outputs[1],
                plots=outputs[0]+'.plots.R',
                num_jobs=options.jobs
            )
    if get_num_files() > 10:
        args += " -an InbreedingCoeff"
        
    run_cmd(gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx16g" % tmp_dir, cpus=options.jobs, mem_per_cpu=16386/options.jobs)


def apply_recalibration_to_snps_or_indels(vcf,recal,tranches,output,mode='SNP',tranche_filter=99.9):
    """Apply the recalibration tranch file to either snps or indels (depending on mode=SNP or mode=INDEL"""
    args = "-T ApplyRecalibration \
            -R {reference} \
            -input {vcf} \
            --ts_filter_level {tranche_filter} \
            -tranchesFile {tranches}  \
            -recalFile {recal} \
            -mode {mode} \
            -nt {num_jobs} \
            -o {output}".format(
                reference=reference,
                vcf=vcf,
                tranche_filter=tranche_filter,
                tranches=tranches,
                recal=recal,
                mode=mode,
                num_jobs=options.jobs,
                output=output)
            
    run_cmd(gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx16g" % tmp_dir)



@follows(find_indel_tranches_for_recalibration)
@files(['multisample.gatk.vcf','multisample.gatk.snp.model','multisample.gatk.snp.model.tranches'],'multisample.gatk.recalibratedSNPS.rawIndels.vcf')
def apply_recalibration_filter_snps(inputs,output):
    apply_recalibration_to_snps_or_indels(inputs[0],inputs[1],inputs[2],output,mode='SNP',tranche_filter=99.9)


@follows(apply_recalibration_filter_snps)
@files(['multisample.gatk.recalibratedSNPS.rawIndels.vcf','multisample.gatk.indel.model','multisample.gatk.indel.model.tranches'],'multisample.gatk.recalibrated.vcf')
def apply_recalibration_filter_indels(inputs,output):
    apply_recalibration_to_snps_or_indels(inputs[0],inputs[1],inputs[2],output,mode='INDEL',tranche_filter=99.0)
    

@transform(apply_recalibration_filter_indels, suffix('.recalibrated.vcf'), '.analysisReady.exome.vcf')
def final_calls(input_vcf, output_vcf):
    """ Produce the final variant calls in the exome regions """    
    args = "-T SelectVariants \
            -R {reference} \
            --variant {input} \
            -o {output} \
            -env -ef \
            -L {exome}".format(
                reference=reference,
                input=input_vcf,
                output=output_vcf,
                exome=exome)

    # remove(input)
    run_cmd(gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx8g" % tmp_dir, mem_per_cpu=8192)


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
            multisample.gatk.recalibratedSNPS.rawIndels.vcf.idx","",
            run_locally=True)
#            multisample.gatk.analysisReady.vcf \


def archive_results():
    # if optional results_archive was not provided - do nothing
    if results_archive == None: return
    
    run_name = os.path.basename(os.getcwd())
    arch_path = os.path.join(results_archive, run_name)
    if not os.path.exists(arch_path): os.mkdir(arch_path)
    run_cmd("cp */*.gatk.bam %s" % arch_path, "", run_locally=True)
    run_cmd("cp */*.gatk.bam.gene_coverage* %s" % arch_path, "", run_locally=True)
    run_cmd("cp */*.exome.vcf %s" % arch_path, "", run_locally=True)
    run_cmd("cp multisample.gatk.gvcf %s" % os.path.join(results_archive,run_name+".multisample.gatk.gvcf"),
            "", run_locally=True)


@posttask(cleanup_files, archive_results)
@follows('final_calls')
@files(split_snp_parameters)
def split_snps(vcf, output, sample):
    """ Split variants by sample, and use sample-specific statistics to filter: AD and DP"""
    AD_threshold=5
    DP_threshold=8
    # what if there is multiple alt alleles??? is the first one most covered      
    args = "-T SelectVariants \
            -R {ref} \
            --variant {vcf} \
            -sn {sample} \
            -select 'vc.getGenotype(\\\"{sample}\\\").getAD().1 >= {ad_thr} && vc.getGenotype(\\\"{sample}\\\").getDP() >= {dp_thr}' \
            -o {out} \
            ".format(ref=reference,
                     vcf=vcf, 
                     sample=sample, 
                     out=output, 
                     ad_thr=AD_threshold, 
                     dp_thr=DP_threshold)
            
    run_cmd(gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx2g" % tmp_dir, mem_per_cpu=2048)


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Main logic


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if __name__ == '__main__':
    if options.just_print:
        pipeline_printout(sys.stdout, options.target_tasks, options.forced_tasks,
                            gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                            verbose=options.verbose, #verbose_abbreviated_path=0,
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
                            multithread     = options.jobs,
                            logger          = stderr_logger,
                            verbose         = options.verbose,
                            gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                            checksum_level  = 0)
    
        
    drmaa_session.exit()
    
