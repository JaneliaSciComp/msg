#!/bin/bash
#
# This installs MSG and all dependencies on a host server.  It is assumed this server will be used for SGE master and all nodes.
# (This is currently set up to run on Ubuntu Linux 12 64 bit. It may not run correctly on other systems.)
# Arguments: 
# 	Arg1 is bin directory to install all binaries and supporting files (make sure this is already on PATH or provide end user 
# 		with a means to get it on their path.
#	Arg2 is msg install directory.  Where to put msg code.

echo "bin install directory:"
echo $1
echo "msg install location:"
echo $2
echo "PATH:"
echo $PATH

#Get MSG software and installers for dependencies and install
cd $2
git clone https://github.com/JaneliaSciComp/msg.git
chmod 777 -R $2
cd msg
make
cd $2/msg/dependencies/

#Install Ghostscript (for R)
apt-get update
apt-get -y install ghostscript

#Install BWA
tar jxvf bwa-0.5.7.tar.bz2 
cd bwa-0.5.7/
make
cp bwa $1
cp *.pl $1
cd $2/msg/dependencies/

#Install samtools
tar jxvf samtools-0.1.9.tar.bz2 
cd samtools-0.1.9/
make
cp -R samtools $1
cp -R misc $1
cp -R bcftools $1
cd $2/msg/dependencies/

#Install BioPython
tar xzvf biopython-1.53.tar.gz 
cd biopython-1.53/
python setup.py build
#python setup.py test
python setup.py install
cd $2/msg/dependencies/

#Install Pyrex
tar xzvf Pyrex-0.9.9.tar.gz
cd Pyrex-0.9.9/
python setup.py install
cd $2/msg/dependencies/

#Install pysam
tar xzvf pysam-0.1.2.tar.gz 
cd pysam-0.1.2/
python setup.py install
cd $2/msg/dependencies/

#Install Stampy
#wget http://www.well.ox.ac.uk/bioinformatics/Software/Stampy-latest.tgz
tar xzvf Stampy-latest.tgz 
cd stampy-1.0.23/
make
#Run this command manually per instructions :-(
g++ `python2.7-config --ldflags` -pthread -shared build/linux-x86_64-2.7-ucs4/pyx/maptools.o build/linux-x86_64-2.7-ucs4/c/maputils.o build/linux-x86_64-2.7-ucs4/c/alignutils.o build/linux-x86_64-2.7-ucs4/readalign.o build/linux-x86_64-2.7-ucs4/algebras.o build/linux-x86_64-2.7-ucs4/frontend.o -o maptools.so
cp -R Stampy/ $1
cp -R plugins/ $1
cp -R maptools.so $1
cp -R ext/ $1
cp stampy.py $1
cd $2/msg/dependencies/

#Install R
apt-get -y install libreadline-dev
tar xzvf R-2.12.2.tar.gz 
cd R-2.12.2/
./configure -prefix=$1/
make
make install
cd $1/bin
mv * ../
cd $1
rm -rf bin

#Install R Libraries
R CMD INSTALL $2/msg/dependencies/HiddenMarkov_1.3-1.tar.gz 
R CMD INSTALL $2/msg/dependencies/zoo_1.6-2.tar.gz 
R CMD INSTALL $2/msg/dependencies/R.methodsS3_1.2.0.tar.gz 
R CMD INSTALL $2/msg/dependencies/R.oo_1.7.3.tar.gz 

#Test dependencies
cd $2/msg
bash test_dependencies.sh
