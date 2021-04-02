MSG: Multiplexed Shotgun Genotyping
http://genomics.princeton.edu/AndolfattoLab/MSG.html
https://github.com/JaneliaSciComp/msg
http://www.ncbi.nlm.nih.gov/pubmed/21233398

# This branch is deprecated!

**Switch to the dev branch instead.**

The master branch contains older, deprectated code.  Use the dev branch of this repository for the latest improvements. Most of the improvements come from the fork at https://github.com/YourePrettyGood/msg, on their dev branch. We are matching their convention rather than merging the changes back to master.

### Dependencies ###
Run test_dependencies.sh to check for the existence of dependencies

NOTE: Most of these packages can be found in the dependencies directory for you to 
extract and install as needed on your system.

Python (2.6.4)
numpy (recent versions should be fine)
bwa (0.5.7)
samtools (0.1.9-3)  ! Newer versions of Samtools will not work !
biopython-1.53
Pyrex-0.9.9
pysam-0.1.2 (apply fix as described here: http://code.google.com/p/pysam/issues/detail?id=22&can=1&q=dandavison0)
R-3.3.1
R packages (HiddenMarkov 1.3-1, zoo 1.6-2, R.methodsS3 3_1.7.1 and R.oo 1.20.0)
Perl Modules (IO::Uncompress::Gunzip)
    -- A command like this should work on Unix Desktop Systems:
    -- perl -MCPAN -e 'install IO::Uncompress::Gunzip'
mailer (Python package: http://pypi.python.org/pypi/mailer) - optional; for sending email alert when run completes

### Installation instructions ###
git clone git://github.com/JaneliaSciComp/msg.git
cd msg
make

### Toy Example ###
The toy example is useful to check that MSG and dependencies are installed and functioning correctly.
TODO

### Example with Real Data ###
Download the data from NCBI's Sequence Read Archive and MSG config and barcodes file
 - example/get_sample_data.sh (requires wget)

