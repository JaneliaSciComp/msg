#!/usr/bin/env python

"""
A command line wrapper around mapping_functions.create_stats
for cases such as the new style parser where we want to run this 
later.

"""

import sys
from mapping_functions import *
from msglib import *
from cmdline.cmdline import CommandLineApp

class Stats(CommandLineApp):
    def __init__(self):
        CommandLineApp.__init__(self)
        
        op = self.option_parser
        op.set_usage('usage: create_stats -i reads.fq -b barcodes')

        op.add_option('-i', '--raw-data', dest='raw_data_file', type='string', default=None, 
                      help='Raw Illumina data file')

        op.add_option('-b', '--barcodes', dest='barcodes_file', type='string', default=None, 
                      help='Barcodes file')

    def main(self):
        print "create_stats.py"
        print "  reads file:", self.options.raw_data_file
        print "  barcodes file:", self.options.barcodes_file
        create_stats(self.options.raw_data_file, self.options.barcodes_file)

if __name__ == '__main__':
    #import doctest; doctest.testmod(); sys.exit()
    try:
        Stats().run()
    except Exception, e:
        print 'Error in create_stats\n'
        print '%s' % e
        sys.exit(2)
