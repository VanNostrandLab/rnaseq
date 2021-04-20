#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A pipeline for processing RNA-Seq data to detect gene differential expression and alternative splicing.
"""

import collections
import glob
import os
import time
import subprocess
import json
import datetime
import sys

import pysam as pysam
import ruffus
import pandas as pd
from loguru import logger

ignored_args = ['target_tasks', 'jobs', 'use_threads', 'just_print', 'log_file', 'verbose', 'forced_tasks']
parser = ruffus.cmdline.get_argparse(description=__doc__, prog='rnaseq',
                                     ignored_args=ignored_args)
parser.add_argument('MANIFEST', type=str,
                    help='Path to the manifest file that specifies paths for RNA-Seq data.')
parser.add_argument('--outdir', type=str,
                    help="Path to the output directory. Default to the current work "
                         "directory and if the specified path does not exist, it will "
                         "try to create it first.",
                    default=os.getcwd())
parser.add_argument('--genome', type=str,
                    help="Path to STAR reference genome index directory.")
parser.add_argument('--repeat', type=str,
                    help="Path to STAR repeat elements index directory.")
parser.add_argument('--gtf', type=str,
                    help="Path to annotation GTF file.")
parser.add_argument('--strand_direction', type=str,
                    help="Read1 strand direction, either forward or reverse.")
parser.add_argument('--strand_specific', type=int, default=2,
                    help='Perform strand-specific read counting. A single integer value, 0 (unstranded), '
                         '1 (stranded) or 2 (reversely stranded) applied to all input files for running '
                         'featureCounts, default: 1 (stranded).')
parser.add_argument('--track', type=str,
                    help="Name for the UCSC Genome Browser track, default: RNA-Seq",
                    default='RNA-Seq')
parser.add_argument('--track_label', type=str,
                    help="Label for the UCSC Genome Browser track, default: RNA-Seq",
                    default='RNA-Seq')
parser.add_argument('--track_genome', type=str,
                    help="Genome name for the UCSC Genome Browser track, default: RNA-Seq",
                    default='hg19')
parser.add_argument('--job', type=str,
                    help="Name of your job, default: RNA-Seq",
                    default='RNA-Seq')
parser.add_argument('--email', type=str,
                    help='Email address for notifying you the start, end, and abort of you job.')
parser.add_argument('--scheduler', type=str,
                    help='Name of the scheduler on your cluster, '
                         'e.g., PBS (or QSUB) or SBATCH (or SLURM), case insensitive.')
parser.add_argument('--time', type=int,
                    help='Time (in integer hours) for running your job, default: 24.',
                    default=24)
parser.add_argument('--memory', type=int,
                    help='Amount of memory (in GB) for all cores needed for your job, default: 32.',
                    default=32)
parser.add_argument('--cores', type=int,
                    help='Number of CPU cores can be used for your job, default: 8.',
                    default=8)
parser.add_argument('--verbose', action='store_true',
                    help='Print out detailed processing messages.')
parser.add_argument('--verbose_path', type=int, default=2,
                    help='Whether input and output paths are abbreviated.')
parser.add_argument('--debug', action='store_true',
                    help='Invoke debug mode.')
parser.add_argument('--dry_run', action='store_true',
                    help='Print out steps and files involved in each step without actually '
                         'running the pipeline.')
parser.add_argument('--hold_submit', action='store_true',
                    help='Generate the submit script but hold it without submitting to the job scheduler. '
                         'Useful when you want to further review the submit script to make sure everything '
                         'looks good and ready to submit.')
parser.add_argument('--target_tasks', metavar='TASKS', action="append", type=str, default=[],
                    help='Target task(s) of pipeline..')
parser.add_argument('--forced_tasks', metavar='TASKS', action="append", type=str, default=[],
                    help='Task(s) which will be included even if they are up to date.')

START_TIME = time.perf_counter()
options = parser.parse_args()
try:
    with open(options.MANIFEST) as f:
        manifest = json.load(f)
except IOError:
    raise ValueError(f'Manifest {options.MANIFEST} may not exist or not be file.')
outdir = manifest.get('outdir', '')
outdir = os.path.dirname(os.path.abspath(options.MANIFEST)) if not outdir or outdir == '.' else ''
setattr(options, 'outdir', options.outdir or outdir)
if os.path.exists(options.outdir):
    if not os.path.isdir(options.outdir):
        raise ValueError(f'Path to outdir {options.outdir} exists, but it is not a directory.')
else:
    os.mkdir(options.outdir)
os.chdir(options.outdir)

genome = options.genome or manifest.get('genome', '')
if not os.path.isdir(genome):
    raise ValueError(f'Reference genome index {genome} is not a directory or does not exist.')
else:
    setattr(options, 'genome', genome)

repeat = options.repeat or manifest.get('repeat', '')
if not os.path.isdir(repeat):
    raise ValueError(f'Repeat element index {repeat} is not a directory or does not exist.')
else:
    setattr(options, 'repeat', repeat)

gtf = options.gtf or manifest.get('gtf', '')
if not os.path.isfile(gtf):
    raise ValueError(f'Annotation GTF {gtf} is not a file or does not exist.')
else:
    setattr(options, 'gtf', gtf)

directions = {'forward': 'f', 'reverse': 'r'}
strand_direction = options.strand_direction or directions.get(manifest.get('strand_direction', ''), '')
if not strand_direction:
    raise ValueError(f'Strand direction {strand_direction} is not a valid direction.')
else:
    setattr(options, 'strand_direction', strand_direction)

setattr(options, 'strand_specific', options.strand_specific or manifest.get('strand_specific', 0))

samples = manifest.get('samples', [])
if not samples:
    raise KeyError('Manifest file does not contain any samples.')
FASTQS, SIZES = {}, []
for sample in samples:
    try:
        fastq1 = sample.get('fastq', '') or sample.get('fastq1', '')
        if os.path.isfile(fastq1):
            SIZES.append(os.path.getsize(fastq1))
        else:
            raise ValueError(f'FASTQ1 {fastq1} is not a file or does not exist.')
        sample['fastq'] = fastq1
        sample['fastq1'] = fastq1
        if 'name' not in sample:
            raise KeyError(f'Sample {sample} does not have name assigned.')
        if 'group' not in sample:
            raise KeyError(f'Sample {sample} does not have group assigned.')
        if 'fastq2' in sample:
            if sample['fastq2']:
                if os.path.isfile(sample['fastq2']):
                    SIZES.append(os.path.getsize(sample['fastq2']))
                else:
                    raise ValueError(f'FASTQ2 {sample["fastq2"]} is not a file or does not exist.')
            else:
                raise ValueError(f'Key fastq2 for sample {sample[name]} does not have a value.')
        link = os.path.join('fastq_to_bam', f'{sample["name"]}.r1.fastq.gz')
        FASTQS[link] = sample
    except KeyError:
        raise KeyError(f'Sample {sample} does not have fastq1 file provided.')

VERBOSE = options.verbose
setattr(options, 'verbose', 3 if options.dry_run else 0)
setattr(options, 'jobs', options.cores)
setattr(options, 'multiprocess', options.cores)
setattr(options, 'use_threads', False)
setattr(options, 'just_print', options.dry_run)
setattr(options, 'verbose_abbreviated_path', options.verbose_path)
setattr(options, 'exceptions_terminate_immediately', True)
setattr(options, 'one_second_per_job', True)
setattr(options, 'log_exceptions', True)
setattr(options, 'logger', ruffus.black_hole_logger)

VERSION = 1.0
logo = fr"""
             _____   _   _
            |  __ \ | \ | |    /\
            | |__) ||  \| |   /  \    ___   ___   __ _
            |  _  / | . ` |  / /\ \  / __| / _ \ / _` |
            | | \ \ | |\  | / ____ \ \__ \|  __/| (_| |
            |_|  \_\|_| \_|/_/    \_\|___/ \___| \__, |
                                                    | |
                           VERSION {VERSION}              |_|
"""

logger.remove()
out = sys.stderr if options.debug else sys.stdout
logger.add(out, format="<light-green>[{time:YYYY-MM-DD HH:mm:ss}]</light-green> <level>{message}</level>",
           filter=lambda record: record["level"].name == "TRACE",
           level="TRACE")
logger.add(out, format="<level>{message}</level>", filter=lambda record: record["level"].name == "DEBUG")
logger.add(out, format="<light-green>[{time:HH:mm:ss}]</light-green> <level>{message}</level>", level="INFO")


@logger.catch(onerror=lambda _: sys.exit(1))
def cmding(cmd, **kwargs):
    """ Run cmd or raise exception if run fails. """
    def format_cmd(command):
        if isinstance(command, str):
            exe = command.split()[0]
        else:
            command = [str(c) for c in command]
            exe = command[0]
            command = ' '.join([f'\\\n  {c}' if c.startswith('-') or '<' in c or '>' in c else c for c in command])
            command = command.splitlines()
            commands = []
            for c in command:
                if len(c) <= 80:
                    commands.append(c)
                else:
                    items = c.strip().replace(' \\', '').split()
                    commands.append(f'  {items[0]} {items[1]} \\')
                    for item in items[2:]:
                        commands.append(' ' * (len(items[0]) + 3) + item + ' \\')
            command = '\n'.join(commands)
            if command.endswith(' \\'):
                command = command[:-2]
        return exe, command
    
    message, start_time = kwargs.pop('message', ''), time.perf_counter()
    program, cmd = format_cmd(cmd)
    if message:
        logger.info(message)
    logger.debug(cmd)
    process = subprocess.Popen(cmd, universal_newlines=True, shell=True, cwd=options.outdir,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
    stdout, stderr = process.communicate()
    if process.returncode:
        logger.error(f'Failed to run {program} (exit code {process.returncode}):\n{stderr or stdout}')
        return process.returncode
    end_time = time.perf_counter()
    run_time = int(end_time - start_time)
    if message:
        message = message.replace(' ...', f' completed in [{str(datetime.timedelta(seconds=run_time))}].')
        logger.info(message)
    return process.returncode


def estimate_process():
    """Estimate number of processes based on the maximum size of fastq file."""
    
    size = max(SIZES) / (1000 * 1000 * 1000) * 4
    n = int(options.memory / size)
    if n == 0:
        n = 1
    elif n > options.cores:
        n = options.cores
    return n


PROCESS = estimate_process()


@ruffus.follows(ruffus.mkdir('fastq_to_bam'))
@ruffus.originate(list(FASTQS.keys()))
def soft_link(link):
    """ Create soft links for original fastq files. """
    
    def make_link(path, link):
        if path:
            if path == os.path.abspath(link):
                message = "No symbolic link was made for {path}! You are directly working on the original file!"
                logger.warning(message)
            else:
                if not os.path.exists(link):
                    message = f'Soft link fastq: {os.path.basename(path)} ...'
                    cmding(f'ln -s {path} {link}', message=message)
    link1, link2 = link, link.replace('.r1.fastq.gz', '.r2.fastq.gz')
    make_link(FASTQS[link]['fastq1'], link1)
    fastq2 = FASTQS[link].get('fastq2', '')
    if fastq2:
        make_link(fastq2, link2)


@ruffus.jobs_limit(1)
@ruffus.transform(soft_link, ruffus.suffix('.fastq.gz'), '.clean.fastq')
def cut_adapt(fastq, output):
    fastq1, fastq2 = fastq, fastq.replace('.r1.fastq.gz', '.r2.fastq.gz')
    output1, output2 = output, output.replace('.r1.', '.r2.')
    cmd = ['cutadapt', '-O', '5', '--times', '2', '-e', '0.0', '-j', options.cores,
           '-m', '18', '--quality-cutoff', '6', '--match-read-wildcards',
           '-b', 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
           '-b', 'AAAAAAAAAAAAAAAAAAAA',
           '-b', 'TTTTTTTTTTTTTTTTTTTT']
    
    if os.path.isfile(fastq2):
        message = f'Cutting adapters for paired reads {fastq1} and\n{45 * " "}{fastq2} ...'
        args = ['-B', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
                '-B', 'TTTTTTTTTTTTTTTTTTTT',
                '-B', 'AAAAAAAAAAAAAAAAAAAA',
                '-o', output1, '-p', output2, fastq1, fastq2]
    else:
        message = f'Cutting adapters for single read {fastq1} ...'
        args = ['-o', output1,  fastq1]
    cmd.extend(args)
    metrics = fastq.replace('.r1.fastq.gz', '.cutadapt.metrics')
    cmd.extend(['>', metrics])
    cmding(cmd, message=message)


@ruffus.jobs_limit(PROCESS)
@ruffus.transform(cut_adapt, ruffus.suffix('.clean.fastq'), '.clean.sort.fastq')
def sort_fastq(fastq, output):
    message = f'Sort fastq {fastq} ...'
    cmd = ['fastq-sort', '--id', fastq, '>', output]
    cmding(cmd, message=message)
    
    fastq2 = fastq.replace('.r1.clean.fastq', '.r2.clean.fastq')
    if os.path.isfile(fastq2):
        output2 = output.replace('.r1.clean.sort.fastq', '.r2.clean.sort.fastq')
        message = f'Sort fastq {fastq2} ...'
        cmd = ['fastq-sort', '--id', fastq2, '>', output2]
        cmding(cmd, message=message)


@ruffus.jobs_limit(1)
@ruffus.follows(ruffus.mkdir('fastq_to_bam/repeat.elements.map'))
@ruffus.transform(sort_fastq,
                  ruffus.formatter(r'.+/(?P<BASENAME>.*).r1.clean.sort.fastq$'),
                  'fastq_to_bam/repeat.elements.map/{BASENAME[0]}/Unmapped.out.mate1')
def map_to_repeat_elements(fastq, mate1):
    fastq1, fastq2 = fastq, fastq.replace('.r1.', '.r2.')
    prefix = os.path.dirname(mate1)
    if not os.path.isdir(prefix):
        os.mkdir(prefix)
    cmd = ['STAR',
           '--runMode', 'alignReads',
           '--runThreadN', options.cores,
           '--alignEndsType', 'EndToEnd',
           '--genomeDir', options.repeat,
           '--genomeLoad', 'NoSharedMemory',
           '--outBAMcompression', 10,
           '--outFileNamePrefix', f"{prefix}/",
           '--outFilterMultimapNmax', 100,
           '--outFilterMultimapScoreRange', 1,
           '--outFilterScoreMin', 10,
           '--outFilterType', 'BySJout',
           '--outReadsUnmapped', 'Fastx',
           '--outSAMattrRGline', 'ID:foo',
           '--outSAMattributes', 'All',
           '--outSAMmode', 'Full',
           '--outSAMtype', 'BAM', 'Unsorted',
           '--outSAMunmapped', 'None',
           '--outStd', 'Log',
           '--readFilesIn', fastq1]
    if os.path.exists(fastq2):
        cmd.append(fastq2)
        message = f'Map paired reads {fastq1} and\n{28 * " "}{fastq2} to repeat elements ...'
    else:
        message = f'Map single read {fastq1} to repeat elements ...'
    cmding(cmd, message=message)


@ruffus.jobs_limit(PROCESS)
@ruffus.transform(map_to_repeat_elements, ruffus.suffix('.out.mate1'), '.out.mate1.sort.fastq')
def sort_mate(mate1, output):
    message = f'Sort mate1 {mate1} ...'
    cmd = ['fastq-sort', '--id', mate1, '>', output]
    cmding(cmd, message=message)
    
    mate2 = mate1.replace('.out.mate1', '.out.mate2')
    if os.path.isfile(mate2):
        output2 = output.replace('.out.mate1.sort.fastq', '.out.mate2.sort.fastq')
        message = f'Sort mate2 {mate2} ...'
        cmd = ['fastq-sort', '--id', mate2, '>', output2]
        cmding(cmd, message=message)


@ruffus.jobs_limit(1)
@ruffus.follows(ruffus.mkdir('fastq_to_bam/reference.genome.map'))
@ruffus.transform(sort_mate,
                  ruffus.formatter(r'.+/Unmapped.out.mate1.sort.fastq'),
                  'fastq_to_bam/reference.genome.map/{subdir[0][0]}/Aligned.out.bam')
def map_to_reference_genome(mate1, bam):
    prefix = os.path.dirname(bam)
    if not os.path.isdir(prefix):
        os.mkdir(prefix)
    mate2 = mate1.replace('Unmapped.out.mate1.sort.fastq', 'Unmapped.out.mate2.sort.fastq')
    cmd = ['STAR',
           '--runMode', 'alignReads',
           '--runThreadN', options.cores,
           '--alignEndsType', 'EndToEnd',
           '--genomeDir', options.genome,
           '--genomeLoad', 'NoSharedMemory',
           '--outBAMcompression', 10,
           '--outFileNamePrefix', f"{prefix}/",
           '--outFilterMultimapNmax', 10,
           '--outFilterMultimapScoreRange', 1,
           '--outFilterScoreMin', 10,
           '--outFilterType', 'BySJout',
           '--outReadsUnmapped', 'Fastx',
           '--outSAMattrRGline', 'ID:foo',
           '--outSAMattributes', 'All',
           '--outSAMmode', 'Full',
           '--outSAMtype', 'BAM', 'Unsorted',
           '--outSAMunmapped', 'None',
           '--outStd', 'Log',
           '--readFilesIn', mate1]
    if os.path.exists(mate2):
        cmd.append(mate2)
        message = f'Map paired mates {mate1} and\n{28 * " "}{mate2} to reference genome ...'
    else:
        message = f'Map single mate {mate1} to reference genome ...'
    cmding(cmd, message=message)


@ruffus.jobs_limit(1)
@ruffus.transform(map_to_reference_genome,
                  ruffus.formatter(r'.+/Aligned.out.bam'),
                  'fastq_to_bam/{subdir[0][0]}.bam')
def sort_index_bam(bam, sorted_bam):
    primary_bam = bam.replace('.bam', '.primary.bam')
    message = f'Creating primary bam from {bam} ...'
    cmd = ['samtools', 'view', '-@', options.cores, '-F', '0x100', '-b', bam, '-o', primary_bam]
    cmding(cmd, message=message)
    
    message = f'Sorting BAM file {bam} ...'
    cmd = ['samtools', 'sort', '-@', options.cores, '-o', sorted_bam, primary_bam]
    cmding(cmd, message=message)
    
    message = f'Indexing BAM file {bam} ...'
    cmd = ['samtools', 'index', sorted_bam]
    cmding(cmd, message=message)


@ruffus.jobs_limit(PROCESS)
@ruffus.follows(ruffus.mkdir('hub'))
@ruffus.transform(sort_index_bam,
                  ruffus.formatter(r'.+/(?P<BASENAME>.*).bam'),
                  'hub/{BASENAME[0]}.plus.bw')
def make_bigwig_files(bam, bigwig):
    def bam_to_bigwig(bam, scale, strand, bw):
        bg, bg_sort = bw.replace('.bw', '.bg'), bw.replace('.bw', '.sort.bg')
        cmd = f'genomeCoverageBed -ibam {bam} -bg -scale {scale} -strand {strand} -du -split > {bg}'
        cmding(cmd)
        cmd = f'bedSort {bg} {bg_sort}'
        cmding(cmd)
        cmd = f'bedGraphToBigWig {bg_sort} {options.genome}/chrNameLength.txt {bw}'
        cmding(cmd)
        cmding(f'rm {bg}')
    
    message, start_time = f'Make BigWig files for {bam} ...', time.perf_counter()
    logger.info(message)
    pos_bw, neg_bw = bigwig, bigwig.replace('.plus.bw', '.minus.bw')
    with pysam.AlignmentFile(bam, 'rb') as sam:
        total_reads = sam.mapped
    r2 = bam.replace('.bam', '.r2.fastq.gz')
    total_reads = total_reads / 2 if os.path.exists(r2) else total_reads
    try:
        scale = 1000000.0 / total_reads
    except ZeroDivisionError:
        logger.error(f'No reads was found in BAM {bam}')
        ruffus.touch_file(bigwig)
        return
    if options.strand_direction in ('f', 'forward'):
        bam_to_bigwig(bam, scale, '+', pos_bw)
        bam_to_bigwig(bam, -1 * scale, '-', neg_bw)
    else:
        bam_to_bigwig(bam, -1 * scale, '-', pos_bw)
        bam_to_bigwig(bam, scale, '+', neg_bw)
    run_time = int(time.perf_counter() - start_time)
    message = message.replace(' ...', f' completed in [{str(datetime.timedelta(seconds=run_time))}].')
    logger.info(message)


@ruffus.merge(make_bigwig_files, 'hub/hub.txt')
def make_hub_files(inputs, output):
    message, start_time = 'Make hub track file ...', time.perf_counter()
    logger.info(message)
    header = f"""hub {options.track.replace(' ', '_')}
shortLabel {options.track_label}
longLabel {options.track_label}
useOneFile on
email {options.email if options.email else 'fei.yuan@bcm.edu'}

genome {options.track_genome}

track {options.track.replace(' ', '_')}
shortLabel {options.track_label}
longLabel {options.track_label}
type bigWig
superTrack on
"""
    block = """
track {basename}
shortLabel {basename}
longLabel {basename}
type bigWig
visibility full
alwaysZero on
autoScale on
aggregate transparentOverlay
showSubtrackColorOnUi on
parent {track}
container multiWig

    track {name1}
    bigDataUrl {plus}
    shortLabel {basename} Plus strand
    longLabel {basename} Plus strand
    type bigWig
    color 0,100,0
    parent {basename}
    
    track {name2}
    bigDataUrl {minus}
    shortLabel {basename} Minus strand
    longLabel {basename} Minus strand
    type bigWig
    color 0,100,0
    parent {basename}
    """
    
    track = options.track.replace(' ', '_')
    with open(output, 'w') as o:
        o.write(header)
        for bw in inputs:
            plus = os.path.basename(bw)
            name1 = plus.replace('.bw', '').replace('.', '_')
            name2 = name1.replace('plus', 'minus')
            basename = plus.replace('.plus.bw', '')
            minus = plus.replace('.plus.bw', '.minus.bw')
            o.write(block.format(track=track, name1=name1, name2=name2, basename=basename, plus=plus, minus=minus))
    run_time = int(time.perf_counter() - start_time)
    message = message.replace(' ...', f' completed in [{str(datetime.timedelta(seconds=run_time))}].')
    logger.info(message)


@ruffus.jobs_limit(1)
@ruffus.follows(make_hub_files)
@ruffus.follows(ruffus.mkdir('differential_expression'))
@ruffus.merge(sort_index_bam, 'differential_expression/count_matrix.tsv')
def feature_count(inputs, output):
    counts = os.path.join('differential_expression', 'feature.count')
    message = 'Running featureCount ...'
    cmd = ['featureCounts', '-a', options.gtf, '-p', '-s', options.strand_specific,
           '-T', options.cores, '-o', counts]
    design, groups = [], []
    for fastq1, sample in FASTQS.items():
        bam = fastq1.replace('.r1.fastq.gz', '.bam')
        cmd.append(bam)
        design.append((f"{sample['name']}_{sample['group']}", sample['group']))
    cmding(cmd, message=message)
    
    logger.info('Parsing feature counts ...')
    start_time = time.perf_counter()
    count_matrix = pd.read_csv(counts, sep='\t', skiprows=1)
    count_matrix = count_matrix.drop(columns=['Chr', 'Start', 'End', 'Strand', 'Length'])
    count_matrix.columns = [os.path.basename(c)[:-4] if c.endswith('.bam') else c
                            for c in count_matrix.columns]
    
    design_matrix = pd.DataFrame(design, columns=['sample', 'group'])
    count_matrix.to_csv(output, index=False, sep='\t')
    design_matrix.to_csv(output.replace('count_matrix.tsv', 'design_matrix.tsv'), index=False, sep='\t')
    end_time = time.perf_counter()
    run_time = int(end_time - start_time)
    logger.info(f'Parsing feature counts complete in [{str(datetime.timedelta(seconds=run_time))}].')


@ruffus.transform(feature_count,
                  ruffus.suffix('count_matrix.tsv'),
                  'significant_result.csv')
def DESeq2(count_matrix, output):
    message = 'Running DESeq2 to identify differentially expressed genes ...'
    cmd = ['deseq.R', count_matrix, count_matrix.replace('count_matrix.tsv', 'design_matrix.tsv')]
    code = cmding(cmd, message=message)
    if code:
        open(output, 'w').close()
        logger.error(f'Running deseq.R failed, empty file {output} was touched.')


@ruffus.jobs_limit(1)
@ruffus.merge(sort_index_bam, 'alternative_splicing/rmats.summary.txt')
@ruffus.follows(ruffus.mkdir('alternative_splicing'))
@ruffus.follows(DESeq2)
def rMATS(inputs, output):
    bams, read_type = collections.defaultdict(list), ''
    for fastq1, sample in FASTQS.items():
        bams[sample['group']].append(fastq1.replace('.r1.fastq.gz', '.bam'))
        fastq2 = sample.get('fastq2', '')
        read_type = 'paired' if fastq2 else 'single'
    b1b2 = []
    for group, files in bams.items():
        file = os.path.join('alternative_splicing', f'{group}_bam_files.txt')
        b1b2.append(file)
        with open(file, 'w') as o:
            o.write(f'{",".join(files)}')
    tmp = os.path.join('alternative_splicing', 'tmp')
    cmd = ['rmats.py',
           '--b1', b1b2[0],
           '--b2', b1b2[1],
           '--gtf', options.gtf,
           '-t', read_type,
           '--readLength', '100',
           '--variable-read-length',
           '--nthread', options.cores,
           '--tmp', tmp,
           '--od', 'alternative_splicing',
           '--libType', 'fr-unstranded',
           '>', 'alternative_splicing/rmats.log']
    message = 'Running rMATS to identify alternative splicing events ...'
    cmding(cmd, message=message)
    os.rename('alternative_splicing/rmats.log', output)


def cleanup():
    logger.info('Deleting soft links ...')
    cmding('rm fastq_to_bam/*.fastq.gz')
    logger.info('Deleting fastq files ...')
    cmding('rm fastq_to_bam/*.clean.fastq')
    logger.info('Compressing sorted fastq files ...')
    for fastq in glob.iglob('fastq_to_bam/*.clean.sort.fastq'):
        cmding(f'pigz -p 8 {fastq}')
    logger.info('Compressing sorted mate files ...')
    for mate in glob.iglob('fastq_to_bam/repeat.elements.map/*/Unmapped.out.mate1'):
        sample = os.path.basename(os.path.dirname(mate))
        cmding(f'pigz -p 8 -c {mate} > fastq_to_bam/{sample}.mate1.sort.fastq.gz')
    for mate in glob.iglob('fastq_to_bam/repeat.elements.map/*/Unmapped.out.mate2'):
        sample = os.path.basename(os.path.dirname(mate))
        cmding(f'pigz -p 8 -c {mate} > fastq_to_bam/{sample}.mate2.sort.fastq.gz')
    logger.info('Deleting map directories ...')
    cmding('rm -rf fastq_to_bam/repeat.elements.map fastq_to_bam/reference.genome.map')
    logger.info('Deleting rMATS temporary directory ...')
    cmding('rm alternative_splicing/rmats.summary.txt')
    cmding('rm -rf alternative_splicing/tmp')


@ruffus.posttask(cleanup)
@ruffus.follows(ruffus.mkdir('qc'))
@ruffus.merge([DESeq2, rMATS, sort_index_bam], 'qc/rnaseq.qc.html')
def qc(inputs, output):
    logger.info('Moving cutadapt metrics to qc ...')
    for metric in glob.iglob('fastq_to_bam/*.cutadapt.metrics'):
        cmding(f'mv {metric} qc/')
    logger.info('Moving STAR logs to qc ...')
    logs = glob.iglob('fastq_to_bam/*/*/Log.final.out')
    for log in logs:
        _, group, basename, _ = log.split('/')
        cmding(f'mv {log} qc/{basename}.{group}.log')
    logger.info('Moving feature count summary to qc ...')
    cmding(f'mv differential_expression/feature.count.summary qc/')
    
    config = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'multiqc.config.yaml')
    cmd = ['multiqc', '--config', config, '--filename', output, '--force', 'qc']
    cmding(cmd, message='Running MultiQC to generating QC summary ...')


def schedule():
    sbatch = """#!/usr/bin/env bash

#SBATCH -n {cores}                        # Number of cores (-n)
#SBATCH -N 1                        # Ensure that all cores are on one Node (-N)
#SBATCH -t {runtime}                  # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem={memory}G                   # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name={job}          # Short name for the job
"""
    sbatch_email = """
#SBATCH --mail-user={email}
#SBATCH --mail-type=ALL
"""
    pbs = """ #!/usr/bin/env bash

#PBS -l nodes=1:ppn={cores}
#PBS -l walltime={runtime}
#PBS -l vmem={memory}gb
#PBS -j oe
#PBS -N {jobname}
"""
    pbs_email = """
#PBS -M {email}
#PBS -m abe
"""
    code = r"""
export TMPDIR={project}/tmp
export TEMP={project}/tmp
export TMP={project}/tmp

{program} \
    --verbose {debug}\
    --outdir {outdir} \
    --genome {genome} \
    --repeat {repeat} \
    --gtf {gtf} \
    --strand_direction {strand_direction} \
    --strand_specific {strand_specific} \
    --cores {cores} \
    {MANIFEST}
"""
    if options.scheduler.upper() in ('PBS', 'QSUB'):
        runtime, directive, exe, mail = f'{options.time}:00:00', pbs, 'qsub', pbs_email
        project = '/project/vannostrand'
    elif options.scheduler.upper() in ('SLURM', 'SBATCH'):
        days, hours = divmod(options.time, 24)
        runtime, directive, exe, mail = f'{days}-{hours:02}:00', sbatch, 'sbatch', sbatch_email
        project = '/storage/vannostrand'
    else:
        raise ValueError(f'Unsupported scheduler: {options.scheduler}, see help for supported schedulers.')
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    setattr(options, 'runtime', runtime)
    setattr(options, 'project', project)
    if options.debug:
        setattr(options, 'debug', '--debug ')
        setattr(options, 'program', os.path.join(root, 'rnaseq', parser.prog))
    else:
        setattr(options, 'debug', '')
        setattr(options, 'program', os.path.join(root, parser.prog))
    text = [directive, mail, code] if options.email else [directive, code]
    text = ''.join(text).format(**vars(options))
    
    submitter = os.path.join(options.outdir, 'submit.sh')
    with open(submitter, 'w') as o:
        o.write(text)
    
    print(f'Job submit script was saved to:\n    {submitter}')
    if options.hold_submit:
        print(f'Job {options.job} was not submitted yet, submit it after carefully review the submit script using:')
        print(f'    {exe} {submitter}')
    else:
        subprocess.run([exe, submitter])
        print(f'Job {options.job} was successfully submitted with the following settings:')
        data = {'Job name:': options.job, 'Output directory:': options.outdir,
                'Number of cores:': options.cores, 'Job memory:': options.memory,
                'Job runtime:': f'{runtime} (D-HH:MM)'}
        for k, v in data.items():
            print(f'{k:>20} {v}')


@logger.catch()
def main():
    if options.scheduler:
        schedule()
    else:
        keys = ('MANIFEST', 'outdir', 'genome', 'repeat', 'gtf', 'strand_direction', 'strand_specific')
        d, strand_specific = vars(options).copy(), {0: '0 (unstranded)', 1: '1 (stranded)', 2: '2 (reversely stranded)'}
        d['strand_direction'] = 'forward' if d['strand_direction'] == 'f' else 'reverse'
        d['strand_specific'] = strand_specific[d['strand_specific']]
        setting = '\n'.join([f'{k:>20}: {v}' for k, v in d.items() if k in keys])
        logger.debug(logo)
        logger.trace(f'Running rnaseq using the following settings:\n{setting}\n')
        try:
            ruffus.cmdline.run(options)
        except OSError:
            pass
        finally:
            logger.debug('')
            run_time = int(time.perf_counter() - START_TIME)
            logger.info(f'Mission Accomplished!')
            logger.trace(f'Time consumed: {str(datetime.timedelta(seconds=run_time))}.')


if __name__ == '__main__':
    main()
