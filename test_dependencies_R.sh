#!/bin/sh

if ! R -f test_dependencies_R.R >/dev/null 2>/dev/null ; then
	echo "Some R packages are missing"
	exit 1
fi

echo "All required R packages are installed"
exit 0
