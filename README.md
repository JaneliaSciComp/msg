# MSG: Multiplexed Shotgun Genotyping
http://genomics.princeton.edu/AndolfattoLab/MSG.html
https://github.com/JaneliaSciComp/msg
http://www.ncbi.nlm.nih.gov/pubmed/21233398

## Dependencies

1. Python (thoroughly tested with 2.7.3 and 2.7.12)
1. NumPy
1. Pyrex
1. BioPython
1. Pysam (should work with any modern version)
1. Perl (any modern version, so typically 5.10+, which will include IO::Uncompress::Gunzip in the base modules)
1. Switch (a Perl module, available as libswitch-perl on Ubuntu)
1. BWA (should work with modern versions, tested with 0.5.7, 0.7.12-r1044, and 0.7.17)
1. SAMtools (**0.1.9-3 is REQUIRED, newer versions do not work**)
1. R (should work with modern versions, tested with 3.2.3, 3.4.1, and 3.5.1, but not 4)
1. R.methodsS3
1. R.oo
1. HiddenMarkov (an R package)
1. A C compiler (GCC and clang preferred, nothing special required)

### Optional dependencies

1. mailer (available from PyPi)

### Known good dependencies from older Janelia versions of MSG

Python (2.6.4)

numpy (recent versions should be fine)

bwa (0.5.7)

samtools (0.1.9-3)  ! Newer versions of Samtools will not work !

biopython-1.53

Pyrex-0.9.9

pysam-0.1.2 (apply fix as described here: http://code.google.com/p/pysam/issues/detail?id=22&can=1&q=dandavison0)

R packages (HiddenMarkov 1.3-1, zoo 1.6-2, R.methodsS3 1.2.0 and R.oo 1.7.3)

Perl Modules (IO::Uncompress::Gunzip)

 * A command like this should work on Unix Desktop Systems: `perl -MCPAN -e 'install IO::Uncompress::Gunzip'`

mailer (Python package: http://pypi.python.org/pypi/mailer) - optional; for sending email alert when run completes

### Installation instructions
```bash
git clone git://github.com/JaneliaSciComp/msg.git
cd msg
#If you want the latest features and fixes, use the dev branch:
git checkout dev
#Install the packaged samtools, and all necessary R packages:
make
```

When you run `make`, note the installation path for the samtools executable that is printed in the output.  You'll need to set `samtools_path` in your msg.cfg to this path.

The Makefile will attempt to install the necessary R packages, and will check for other major dependencies (Perl, Python, R, BWA, SAMtools, BioPython, Pysam)

### Docker-based installation (non-cluster usage)

Thanks to @ychenbioinfo for doing the heavy lifting getting a Dockerfile set up for MSG!

```bash
docker pull youreprettygood/msg
docker run -i -t --name [Name for MSG container instance] -v [Absolute path to MSG working directory]:/data:Z youreprettygood/msg:0.5.1 /bin/bash
#If you need to use reads to update your parental genomes,
# you can use MSG to do so (although better methods exist)
# by calling:
msgUpdateParentals.pl
#When you have your working directory all set (including parental references),
# run MSG:
msgCluster.pl
#Note that the current Dockerfile does not support cluster job submission,
# so you MUST have cluster=0 in update.cfg and msg.cfg
```

## Toy Example

The toy example is useful to check that MSG and dependencies are installed and functioning correctly.

To test your installation, run MSG (according to the instructions above, e.g. the Docker instructions) in the `example_MSG_toy` directory, and then compare the `ancestry-prob-par*.tsv` files to those found in `example_results` under the appropriate version of BWA.

Be sure to adjust `update.cfg` and `msg.cfg` if you wish to test a cluster installation. It's best to run `msgUpdateParentals.pl` on a single node (i.e. `cluster=0`), but you can change `msg.cfg` to use `cluster=1` as long as you adjust other cluster job submission options appropriately.

The toy example requires that you run `msgUpdateParentals.pl` before you run `msgCluster.pl`.

Once `msgUpdateParentals.pl` and `msgCluster.pl` have successfully completed, there are a few different primary output files to check: `ancestry-prob-par*.tsv` (the full matrices of posterior probabilities of ancestry assignments), the per-individual HMM plots `hmm_fit/indiv*/indiv*-hmmprob.pdf` (showing the posterior probability of each ancestry along each scaffold), and the diagnostic plots in `hmm_fit_images/`.

The primary test for a successful installation is that there are no differences between the `ancestry-prob-par*.tsv` files from your run and the versions found in `example_results` under your version of BWA. We provide a simple script `diff_posteriors.awk` to identify any differences between the corresponding files. You can check this with:

```bash
[path to MSG install]/diff_posteriors.awk ancestry-prob-par1.tsv example_results/bwa_[version]/ancestry-prob-par1.tsv
[path to MSG install]/diff_posteriors.awk ancestry-prob-par1par2.tsv example_results/bwa_[version]/ancestry-prob-par1par2.tsv
[path to MSG install]/diff_posteriors.awk ancestry-prob-par2.tsv example_results/bwa_[version]/ancestry-prob-par2.tsv
```

Thus, a successful installation and test will have zero output from each of these three commands.

You can visually inspect the per-individual HMM plots and diagnostic plots if you want, though identical `ancestry-prob-par*.tsv` files should imply identical per-individual HMM plots.

## Example with Real Data

**Instructions incomplete, pending update**

Download the data from NCBI's Sequence Read Archive and MSG config and barcodes file

 * example/get_sample_data.sh (requires wget)

## Potential future example datasets (not included yet, slow and take lots of disk space):

1. [Cande et al. (2012) PLoS One](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0043888)

2. [Liu et al. (2019) Current Biology](https://www.sciencedirect.com/science/article/abs/pii/S0960982219306888)

## Options

### Dependency resolution:

* `samtools_path`: (path, Default=`samtools`) Path to the samtools 0.1.9 binary. Default is to use the samtools version in your PATH, but beware that MSG requires samtools 0.1.9 and will fail with newer versions

### Input files:

* `parent1` and `parent2`: (path, no default) Path to parental (par1 and par2) reference FASTAs
* `parent1_reads` and `parent2_reads`: (path, no default) Path to parental (par1 and par2) DNAseq reads to use for updating variant sites in the references
* `barcodes`: (path, no default) Path to barcodes file (TSV, 4 columns: barcode sequence, sample ID, plate ID, sex). Sex is interpreted as if the species was male-heterogametic, so you may need to flip your sex labels if your species is female-heterogametic
* `reads`: (path, no default) Path to reads FASTQ file. If you already parsed your per-sample reads, create an empty file with this name, and place your per-sample FASTQs into a subdirectory named this (plus `_parsed`). The per-sample FASTQs must be named following this format: `indiv[sampleID]_[sampleBarcode]`. If the filename ends in `.gz`, MSG will assume the reads are gzipped

### Parsing:

* `re_cutter`: (string, Default=MseI) Restriction enzyme used for library preparation (MseI was used in the original Andolfatto et al. 2011 paper). Ignored if the `${reads}_parsed` directory already exists
* `linker_system`: (string, Default=Dros_SR_vII) Linker system used for restriction-based libraries. This determines how adapters are extracted and parsed from the input reads
* `new_parser`: (0 or 1, Default=0) Use new read parser?
* `new_parser_offset`: ?
* `new_parser_filter_out_seq`: ?
* `index_file`: (path, Default='') 
* `index_barcodes`: (path, Default='') 

### Read quality trimming:

* `quality_trim_reads_thresh`: (non-negative integer, Default=0) Threshold PHRED-scaled base quality for trimming
* `quality_trim_reads_consec`: (non-negative integer, Default=30) Minimum number of bases with base qualities below `quality_trim_reads_thresh` required for trimming

### Read mapping:

* `bwaindex1` and `bwaindex2`: (bwtsw, is, rb2, or auto, Default=bwtsw) BWA index type to build for each parental reference.  Typically bwtsw
* `bwa_alg`: (aln, bwasw, or mem, Default=aln) BWA mapping method to use, aln and bwasw are well-tested, mem is experimental
* `bwa_threads`: (positive integer, Default=1) Number of threads to use for BWA mapping. Typically 8
* `use_stampy`: (0 or 1, Default=0) Whether or not to use STAMPY for mapping
* `indiv_stampy_substitution_rate`: (float, default=0.001) Substitution rate for STAMPY
* `stampy_premap_w_bwa`: (0 or 1, Default=1) Whether or not to use a BWA step to speed up STAMPY. Only relevant when `use_stampy = 1`
* `stampy_pseudo_threads`: (positive integer, Default=0) Number of pseudo-threads STAMPY will use for mapping. Only applies when `cluster = 1`

### Read and AIM filters:

* `max_mapped_reads`: (non-negative integer, Default=0) Maximum number of reads per individual to use, applied by truncating SAM. Setting to 0 means no maximum
* `one_site_per_contig`: (0 or 1, Default=1) Use only one AIM per region of contiguous non-zero depth (plus fuzz of 5 bp) within an individual. This is supposed to mitigate issues of correlated evidence from multiple AIMs in a single read
* `filter_hmm_data`: (0 or 1, Default=0) Filter AIMs using `filter_hmmdata.pl` instead of via `one_site_per_contig` (be sure to set `one_site_per_contig=0`). This filter is less stringent for higher depth datasets, since regions of contiguous non-zero depth are more frequent even though evidence for proximal AIMs may be derived from independent reads
* `read_length`: (positive integer, Default=100) Typical read length to use for the `filter_hmm_data` heuristic. This heuristic is fairly simple: Eliminate the second AIM of any pair separated by less than `read_length`
* `AS_XS_threshold`: (positive integer, Default=6) Difference between primary and secondary alignment scores from BWA-MEM for a read to be considered uniquely mapping and retained
* `thinfac`: (float > 0.0 and <= 1.0, Default=1.0) Fraction of AIMs to keep, sampled at random. Applied during msgRun3 for the ancestry-prob-par*-thinned.tsv output file
* `difffac`: (float > 0.0 and <= 1.0, Default=0.01) Euclidean difference between vectors of posterior probabilities of adjacent AIMs
* `pnathresh`: (float > 0.0 and <= 1.0, Default=0.03) Maximum fraction of missing data at an AIM to retain when computing the LOD matrix
* `pepthresh`: ?
* `gff_thresh_conf`: ?
* `theta`: (float > 0.0 and <= 1.0, Default=1.0) Scaling factor to account for dependence between reads

### HMM tuning:

* `priors`: (3 comma-separated floats, Default=0,0.5,0.5) Prior probabilities of the three ancestry states: Homozygous par1 (par1/par1), Heterozygous (par1/par2), and Homozygous par2 (par2/par2), respectively. Must sum to 1.0. Default is for an F1 backcross to par2. An F2 intercross would have `priors=0.25,0.5,0.25`
* `deltapar1` and `deltapar2`: (float, Default=0.01) Error rates for bases in the parental (par1 and par2, respectively) references. Typically these are tuned to match the mode of the error_gamma plots
* `recRate`: (float, Default=0) Total number of expected recombination events per generation across the entire genome. Typically a good starting point is the number of chromosomes in your species (e.g. 3 for many Drosophila species)
* `rfac`: (float > 0.0 and <= 1.0, Default=0.00001) Scaling factor for `recRate` to reduce the overall HMM switching rate

### Chromosome/scaffold handling:

* `chroms`: (comma-separated strings, Default=all) Chromosomes to process through the HMM. Default of `all` means to use all scaffolds in the par2 reference
* `sexchroms`: (comma-separated strings, Default=X) Chromosome(s) to treat as hemigametic in samples labeled as "male". Haplodiploid species can specify `all` here
* `chroms2plot`: (comma-separated strings, Default=all) Chromosomes to display in the HMM plots. Typically this should be the longest scaffolds to avoid messy plots due to including many short scaffolds

### Output options:

* `full_summary_plots`: (0 or 1, Default=1) If 1, run `summaryPlots.R` to generate full diagnostic plots. If 0, run `combine.py` (may have bugs)
* `plot_lod_matrix`: (0 or 1, Default=1) Whether or not to generate the LOD matrix, which is the most memory-intensive step of msgRun3. Memory usage is O(m^2) where m is the number of AIMs across all individuals, which increases approximately logarithmically with the number of individuals. This can reach absurd amounts of RAM (~100-200 GB) for about a thousand Drosophila individuals

### Parental reference updating:

* `update_minQV`: (non-negative integer, no default) Minimum consensus quality PHRED score required for updating a reference base with the consensus base in the pileup
* `min_coverage`: (positive integer, no default) Minimum number of reads covering a site to allow updating
* `max_coverage_stds`: (float, no default) Maximum number of standard deviations from the mean that depth is allowed to vary for a site to be updated
* `max_coverage_exceeded_state`: (N or D, no default) Rule used for updating sites exceeding maximum coverage. N means mask, D means default(?)

### Partial pipeline execution:

* `msgRun1`: (0 or 1, Default=1) Whether to run msgRun1 (preparation of parental references for mapping)
* `msgRun2`: (0 or 1, Default=1) Whether to run msgRun2 (mapping, filtering, fitting of HMM for each individual)
* `msgRunOther`: (0 or 1, Default=1) Whether to run msgRun2a, msgRun3, msgRun4, msgRun5, msgRun6 (msgRun3 merges and interpolates posterior probabilities across individuals to make the ancestry-prob-par*.tsv files, and some diagnostics; the other jobs do other minor tasks)

### Parallelization and cluster submission:

* `n_parallel`: (positive integer, Default=1) Degree of parallelization of parent1or2-hmm.sh calls (second half of msgRun2). This is implemented using the fork(2) system call and a named pipe (FIFO) acting as a semaphore
* `cluster`: (0 or 1, Default=1) Whether or not MSG should be run by submitting jobs to a cluster
* `threads`: (positive integer, Default=8) 
* `barcodes_per_job`: (positive integer, Default=1) Number of barcodes per msgRun2 task within the cluster job array. Setting to 1 will likely result in the fastest run if the cluster queue is empty and there are many available compute nodes, 10 is more sane and behaves better when the queue is packed

SGE-related:

* `addl_qsub_option_for_exclusive_node`: (string, Default='') 
* `addl_qsub_option_for_pe`: (string, Default='') 
* `custom_qsub_options_for_all_cmds`: (string, Default='') 

Generic to clusters:

* `submit_cmd`: (string, Default=`qsub -N $jobname -cwd $params{"addl_qsub_options_for_exclusive_node"}$params{"custom_qsub_options_for_all_cmds"}-b y -V -sync n`) Cluster job submission command with base options (e.g. queue/partition, QOS, etc.). Default is for SGE, SLURM would use `sbatc -J ${jobname} -o ${logdir}/${jobname}.%j.stdout -e ${logdir}/${jobname}.%j.stderr`, LSF would use `bsub -J ${jobname} -o ${logdir}/${jobname}.stdout -e ${logdir}/${jobname}.stderr`
* `default_submit_options`: (string, Default='') Options for submission of msgRun1, 2a, 2b, and 4-7
* `msgrun2_submit_options`: (string, Default='') Options for submission of msgRun2
* `msgrun3_submit_options`: (string, Default='') Options for submission of msgRun3
* `depend_arg`: (string, Default=`-hold_jid $prev_jobname`) Argument to set dependency of current (non-array) job on previous job. Default is for SGE, SLURM would use `-d afterok:${prev_jobid}`, LSF would use `-w "done(${prev_jobid})"`
* `array_job_variable`: (string, Default=`$SGE_TASK_ID`) Environment variable name for the current task ID within the job array. Default is for SGE, SLURM would use `$SLURM_ARRAY_TASK_ID`, LSF would use `$LSB_JOBINDEX`
* `array_job_arg`: (string, Default=`-t 1-${array_size}`) Argument to submit a job array. Default is for SGE, SLURM would use `-a 1-${array_size}`, LSF would use `-J "${jobname}[1-${array_size}]"`
* `array_job_depend_arg`: (string, Default=`-hold_jid $prev_jobname`) Argument to set dependency of current array job on previous job. Default is for SGE, SLURM would use ``, LSF would use ``

### Email notifications:

* `notify_emails`: (string, Default='')
* `email_host`: (string, Default='') 

### Verbosity/logging:

* `debug`: (0 or 1, Default=0) Extra logging
* `verbose`: (0 or 1, Default=0) Extra logging

## TODO:

[] Generate a Singularity image for deployment on a cluster
[] Generalize cluster submission in `msgUpdateParentals.pl`
[] Figure out a way to reduce the memory usage when generating the LOD matrix
[] Eliminate (if possible) the samtools version dependency
