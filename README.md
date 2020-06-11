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
docker run -i -t --name [Name for MSG container instance] -v [Absolute path to MSG working directory]:/data:Z youreprettygood/msg:0.5 /bin/bash
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

## Potential future example datasets

1. [Cande et al. (2012) PLoS One](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0043888)

2. [Liu et al. (2019) Current Biology](https://www.sciencedirect.com/science/article/abs/pii/S0960982219306888)
