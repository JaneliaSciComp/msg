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

if ! which samtools >/dev/null 2>/dev/null ; then
	echo "Please install samtools somewhere in your PATH"
	exit 1
fi

echo "All required executables found in PATH"
exit 0
