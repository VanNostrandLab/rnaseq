# rnaseq

A pipeline for processing RNA-Seq data to detect gene differential expression and alternative splicing.

## Usage:
On taco cluster, issue the following command to check the usage:

```shell script
$ /storage/vannostrand/software/rnaseq/rnaseq -h
usage: rnaseq [-h] [--version] [--touch_files_only] [--recreate_database]
              [--checksum_file_name FILE] [--flowchart FILE]
              [--key_legend_in_graph] [--draw_graph_horizontally]
              [--flowchart_format FORMAT] [--outdir OUTDIR] [--genome GENOME]
              [--repeat REPEAT] [--gtf GTF]
              [--strand_direction STRAND_DIRECTION]
              [--strand_specific STRAND_SPECIFIC] [--job JOB] [--email EMAIL]
              [--scheduler SCHEDULER] [--time TIME] [--memory MEMORY]
              [--cores CORES] [--verbose] [--dry_run] [--hold_submit]
              [--target_tasks TASKS] [--forced_tasks TASKS]
              MANIFEST

A pipeline for processing RNA-Seq data to detect gene differential expression
and alternative splicing.

positional arguments:
  MANIFEST              Path to the manifest file that specifies paths for
                        RNA-Seq data.

optional arguments:
  -h, --help            show this help message and exit
  --outdir OUTDIR       Path to the output directory. Default to the current
                        work directory and if the specified path does not
                        exist, it will try to create it first.
  --genome GENOME       Path to STAR reference genome index directory.
  --repeat REPEAT       Path to STAR repeat elements index directory.
  --gtf GTF             Path to annotation GTF file.
  --strand_direction STRAND_DIRECTION
                        Read1 strand direction, either forward or reverse.
  --strand_specific STRAND_SPECIFIC
                        Perform strand-specific read counting. A single
                        integer value, 0 (unstranded), 1 (stranded) or 2
                        (reversely stranded) applied to all input files for
                        running featureCounts, default: 1 (stranded).
  --job JOB             Name of your job, default: chimeras
  --email EMAIL         Email address for notifying you the start, end, and
                        abort of you job.
  --scheduler SCHEDULER
                        Name of the scheduler on your cluster, e.g., PBS (or
                        QSUB) or SBATCH (or SLURM), case insensitive.
  --time TIME           Time (in integer hours) for running your job, default:
                        24.
  --memory MEMORY       Amount of memory (in GB) for all cores needed for your
                        job, default: 32.
  --cores CORES         Number of CPU cores can be used for your job, default:
                        8.
  --verbose             Print out detailed processing messages.
  --dry_run             Print out steps and files involved in each step
                        without actually running the pipeline.
  --hold_submit         Generate the submit script but hold it without
                        submitting to the job scheduler. Useful when you want
                        to further review the submit script to make sure
                        everything looks good and ready to submit.
  --target_tasks TASKS  Target task(s) of pipeline..
  --forced_tasks TASKS  Task(s) which will be included even if they are up to
                        date.

Common options:
  --version             show program's version number and exit

pipeline arguments:
  --touch_files_only    Don't actually run any commands; just 'touch' the
                        output for each task to make them appear up to date.
  --recreate_database   Don't actually run any commands; just recreate the
                        checksum database.
  --checksum_file_name FILE
                        Path of the checksum file.
  --flowchart FILE      Don't run any commands; just print pipeline as a
                        flowchart.
  --key_legend_in_graph
                        Print out legend and key for dependency graph.
  --draw_graph_horizontally
                        Draw horizontal dependency graph.
  --flowchart_format FORMAT
                        format of dependency graph file. Can be 'svg', 'svgz',
                        'png', 'jpg', 'psd', 'tif', 'eps', 'pdf', or 'dot'.
                        Defaults to the file name extension of --flowchart
                        FILE.
```

Note: when you are not inside the rnaseq directory, call the executable script `rnaseq` using 
its relative path or absolute path.