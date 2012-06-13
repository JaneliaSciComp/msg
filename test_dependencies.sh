#!/bin/sh

sh test_dependencies_executable.sh || exit 1
sh test_dependencies_R.sh || exit 1
python test_dependencies_python.py || exit 1
perl test_dependencies_perl.pl || exit 1
exit 0
