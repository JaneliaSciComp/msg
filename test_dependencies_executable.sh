#!/bin/sh

if ! perl --version >/dev/null 2>/dev/null ; then
	echo "Please install perl somewhere in your PATH"
	exit 1
fi

if ! python -V >/dev/null 2>/dev/null ; then
	echo "Please install python somewhere in your PATH"
	exit 1
fi

if ! R --version >/dev/null 2>/dev/null ; then
	echo "Please install R somewhere in your PATH"
	exit 1
fi

if ! which bwa >/dev/null 2>/dev/null ; then
	echo "Please install bwa somewhere in your PATH"
	exit 1
fi

samtools_required_version="Version: 0.1.9 (r783)"
if ! [ -e $(dirname $0)/samtools ] ; then
        echo "Please install samtools ($samtools_required_version) within the msg directory"
	    exit 1
fi
samtools_version=`$(dirname $0)/samtools 2>&1 | grep 'Version:'`
if [ "$samtools_version" != "$samtools_required_version" ] ; then
    echo "The samtools found within the msg directory is $samtools_version.  It must be $samtools_required_version" 
    exit 1
fi

echo "All required executables found in PATH"
exit 0
