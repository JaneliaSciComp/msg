#!/bin/sh

if ! [ -x "$(command -v perl)" ]; then
#if ! perl --version >/dev/null 2>/dev/null ; then
	echo "Please install perl somewhere in your PATH"
	exit 1
fi

if ! [ -x "$(command -v python)" ]; then
#if ! python -V >/dev/null 2>/dev/null ; then
	echo "Please install python somewhere in your PATH"
	exit 1
fi

if ! [ -x "$(command -v Rscript)" ]; then
#if ! R --version >/dev/null 2>/dev/null ; then
	echo "Please install R somewhere in your PATH"
	exit 1
fi

if ! [ -x "$(command -v bwa)" ]; then
#if ! which bwa >/dev/null 2>/dev/null ; then
	echo "Please install bwa somewhere in your PATH"
	exit 1
fi

if (! [ -x "$(command -v samtools)" ]) && (! [ -x "samtools-0.1.9/samtools" ]); then
	echo "Please install samtools 0.1.9 that came with MSG using the Makefile"
	exit 1
fi

echo "All required executables found in PATH"
exit 0
