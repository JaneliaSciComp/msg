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
1. BWA (should work with modern versions, tested with 0.5.7 and 0.7.12-r1044)
1. SAMtools (**0.1.9-3 is REQUIRED, newer versions do not work**)
1. R (should work with modern versions, tested with 3.2.3 and 3.4.1)
1. R.methodsS3
1. R.oo
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
    -- A command like this should work on Unix Desktop Systems:
    -- perl -MCPAN -e 'install IO::Uncompress::Gunzip'
mailer (Python package: http://pypi.python.org/pypi/mailer) - optional; for sending email alert when run completes

### Installation instructions ###
```bash
git clone git://github.com/JaneliaSciComp/msg.git
cd msg
make
```

When you run `make`, note the installation path for the samtools executable that is printed in the output.  You'll need to set `samtools_path` in your msg.cfg to this path.

The Makefile will attempt to install the necessary R packages, and will check for other major dependencies (Perl, Python, R, BWA, SAMtools, BioPython, Pysam)

## Toy Example
The toy example is useful to check that MSG and dependencies are installed and functioning correctly.
**Instructions incomplete, pending update**

## Example with Real Data
**Instructions incomplete, pending update**
Download the data from NCBI's Sequence Read Archive and MSG config and barcodes file
 - example/get_sample_data.sh (requires wget)

