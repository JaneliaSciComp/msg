#!/usr/bin/env python 
from datetime import datetime
from optparse import OptionParser
import fileinput
import gzip
import re
import sys

__doc__ = """
TQS

Trim Quality Sequences (TQS)

SYNOPSIS
   Quality trim FASTQ sequence reads using user-defined thresholds 
"""
__author__ = "Rene L. Warren"
__version__ = '1.2'

#LICENSE
#   Copyright (c) 2007 Canada's Michael Smith Genome Science Centre.  All rights reserved.

#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

# Modified by Lance Parsons at Princton University's Lewis-Sigler Institute for Integrative Genomics
#    Adapted to trim "standard" FASTQ files (PHRED+33) (6/18/2009)
#    Adapted to output FASTA or FASTQ files (1/25/2010)
#    Reduce memory footprint, do not load entire file, no variable grows (1/25/2010)
#    Fixed bug when seq_len < consecutive (1/25/2010)
#    Allow reading of gzipped files directly (8/18/2010)
#    Allow output of gzipped files (12/10/2010)
#    Allow specification of output name (12/17/2010)
# Greg Pinero
#   Output file is output_base + .trim.fast(a/q) (.gz) (2/23/2011)



def main():
    usage = "Usage: %s --help"

    parser = OptionParser()
    parser.add_option("-f", "--fastq file", dest="fastqfile", help="Sanger encoded fastq file - PHRED quality scores, ASCII+33",)
    parser.add_option("-t", "--Phred quality threshold", dest="threshold", type="int", default=10, help="Base intensity threshold value (Phred quality scores 0 to 40, default=10)",)
    parser.add_option("-c", "--consec", dest="consec", type="int", default=20, help="Minimum number of consecutive bases passing threshold values (default=20)",)
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", help="Runs in Verbose mode.",)
    parser.add_option("-q", "--qualities", dest="qualities", action="store_true", help="Outputs Qualities to FASTQ file (default is FASTA)",)
    parser.add_option("-z", "--zip", dest="gzip", action="store_true", default=False, help="Compress output with gzip",)
    parser.add_option("-o", "--output", dest="output_base", default="", help="Output filename base",)
    (opts, args) = parser.parse_args()
    if len(args) != 0:
        parser.error("incorrect number of arguments")

    try:
        #f = open(opts.fastqfile)
        f = fileinput.FileInput(opts.fastqfile, openhook=fileinput.hook_compressed)
        # seq = f.readlines()
        # f.close()
    except Exception, e:
        print "ERROR: Could not read from %s: %s" % (opts.fastqfile, e)
        print usage % (sys.argv[0:])
        sys.exit()

    output_base = opts.fastqfile
    if opts.output_base != "":
        output_base = opts.output_base
    if opts.qualities:
        fasta = "%s.trim.fastq" % (output_base)
    else:
        fasta = "%s.trim.fasta" % (output_base)
    
    if opts.gzip:
        fasta = "%s.gz" % fasta
        
    log = "%s_T%sC%s.trim.log" % (output_base, opts.threshold, opts.consec)
    minimum_length = 15
    
    try:
        if opts.gzip:
            FASTA = gzip.open(fasta, 'wb', 9)
        else:
            FASTA = open(fasta, 'w')
    except:
        print "ERROR: Can not write to %s" % fasta
        sys.exit()

    try:
        LOG = open(log, 'w')
    except:
        print "ERROR: Can not write to %s" % log
        sys.exit()
    
    if opts.consec < minimum_length:
        print "ERROR: -c must be a number larger than %i." % (minimum_length)
        sys.exit()
        
    #print h.heap()

    LOG.write("""
Running:
%s
-f %s
-c %s
-t %s
Fasta file: %s

""" % (sys.argv[0:], opts.fastqfile, opts.consec, opts.threshold, fasta))
    
    t1 = datetime.now()
    LOG.write("\n\nTrimming low quality bases: %s\n" % str(t1)[:len('2006-10-05 23:04')])
    readNtrim(f, opts.threshold, opts.consec, opts.verbose, opts.qualities, FASTA, LOG)
    LOG.write("DNA sequences have been trimmed accordingly and placed in %s" % fasta)
    
    #print h.heap()
    f.close()
    LOG.close()
    FASTA.close()
    
    #print h.heap()
    return    

#--------------------------------------------------------------------------------------
def readNtrim(fastq, threshold, consecutive, verbose, qualities, FASTA, LOG):
    """
    Parse a solexa-illumina export file
    SOLEXA3_77_30V9CAAXX        4    1    1068    522        1    GGACAGCTGACAGCTGTTAAGAAGGACCCTATGTTAAAGGAAATGGATAC    YYYYYYYYYYYJYY
YYYYRYYYYYYYYYYYTTTTTOOOMOOOMMOOOOOG    chr13        36311743    F    50    52    121            187    R    N
    Return a Dictionary of sequence order number, with the index value and length to extract 
    """
    #trim_info = {}
    ok_read = 0
    read_number = 0
    record_line = 0
    
    if verbose:
        print "Printing trimming pattern for all reads passing the set threshold values...\n"
    
    for line in fastq:
        record_line += 1
        if record_line == 1:
            read_id = line.strip()
        elif record_line == 2:
            seq = line.strip()
        elif record_line == 3:
            qual_id = line.strip()
            assert (read_id[1:] == qual_id[1:]) or (qual_id[1:] == '')
        elif record_line == 4:
            record_line = 0
            qual = line.strip()
            read_number += 1
            concat = ""            ### concat builds a string of bases passing the user-defined filter 
            """
            print "line%s\tseq:%s\tqual:%s\n" % (line,info[6],info[7])
            """
            pos = 0
            for qual_char in qual:
                Q = (ord(qual_char) - 33)
                pos += 1
                if Q < threshold:
                    concat += "x"
                else:
                    concat += "-"
                """
                print "base#%i. Illumina qual (%s) == phredQ (%i)\n" % (pos,illumina_qual,Q)
                """
    
            seq_len = len(seq)
            head_match_regex = re.compile(r"\-{%i,%i}" % (consecutive, max(consecutive, seq_len))) 
            head_match = head_match_regex.search(concat)
            if head_match != None:
                ok_read += 1
                col = head_match.span()
                            #if not trim_info.has_key(read_number):
                            #        trim_info[read_number] = {}
    
                start = int(col[0])    
                end = int(col[1])
                
                trim_seq = seq[start:end]
                if qualities:
                    trim_qual = qual[start:end]
                    FASTA.write("%s\n%s\n+\n%s\n" % (read_id, trim_seq, trim_qual))
                else:
                    FASTA.write(">%s\n%s\n" % (read_id, trim_seq))
    
                if verbose:
                    print "passed seqs:%i line#%i %s (start trim:%i,end trim:%i) %s" % (ok_read, read_number, concat, start, end, trim_seq)
            
            if read_number % 100000 == 0:
                #print h.heap()
                print "Read %i sequences" % (read_number)
            
    LOG.write("%i out of %i sequences passed your filter (-t >= %i and -c >= %i)\n" % (ok_read, read_number, threshold, consecutive))

    return



if __name__ == '__main__':
    main()
    sys.exit()
