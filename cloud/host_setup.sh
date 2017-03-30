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

#Install s3cmd (for accessing Amazon S3 hosted files)
cd /tmp
wget https://github.com/s3tools/s3cmd/archive/v1.5.0-rc1.tar.gz
tar xzvf v1.5.0-rc1.tar.gz 
cd s3cmd-1.5.0-rc1/
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

# Install R Manually for now (Or contact Patrick for a better install script?)
# This approach worked on Ubuntu Precise
#codename=$(lsb_release -c -s)
#echo "deb http://cran.stat.ucla.edu/bin/linux/ubuntu $codename/" | sudo tee -a /etc/apt/sources.list > /dev/null
#sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
#sudo add-apt-repository ppa:marutter/rdev
#sudo apt-get update
#sudo apt-get upgrade
#sudo apt-get install r-base r-base-dev
##N for Configuration file `/etc/cloud/cloud.cfg' question	
#sudo apt-get install r-base-core

#Install R Libraries
# (March 2017) This script worked well to install (zoo no longer required)
# https://raw.githubusercontent.com/YourePrettyGood/msg/master/check_dependencies_R.R

#Test dependencies (Won't pass without R)
cd $2/msg
bash test_dependencies.sh
