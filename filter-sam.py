#!/usr/bin/env python

"""
This script removes any reads from the SAM file that were not mapped.
For BWA it also keeps only reads with one best match and no suboptimal matches.

"""

import sys
from cmdline.cmdline import CommandLineApp
import pysam

class SamFilter(CommandLineApp):
    def __init__(self):
        CommandLineApp.__init__(self)
        op = self.option_parser
        op.set_usage('usage: filter-sam.py')
        op.add_option('-i', '--infile', dest='infile', type='string', default=None, 
                      help='Input .sam file')
        op.add_option('-o', '--outfile', dest='outfile', type='string', default=None, 
                      help='Output .sam file')
        op.add_option('-a', '--bwa_alg', dest='bwa_alg', type='string', default='', 
                      help='Output .sam file')
        op.add_option('-s', '--use_stampy', dest='use_stampy', type='int', default=0, 
                      help='SAM file was mapped with stampy')

    def main(self):
        #my best guess about what the SAM flags mean:
        flags_to_consider = set([
            0,  #no information
            4,  #segment unmapped
            16, #SEQ being reverse complemented
        ])
    
        bwa_alg = self.options.bwa_alg
        use_stampy = self.options.use_stampy
        infile = pysam.Samfile(self.options.infile, 'r')
        outfile = pysam.Samfile(self.options.outfile, 'wh', template=infile)
        #outfile = pysam.Samfile(self.options.outfile, 'w', template=infile)
        print('Filtering reads in %s to %s' % (infile, outfile))
        i = 0
        omit = 0
        unmapped = 0
        for read in infile.fetch():
            i += 1
  
            if not (read.flag in flags_to_consider):
                print 'read %s flag %d not in {0,4,16}: omitting' % (read.qname, read.flag)
                continue
            #skip unmapped reads
            if read.flag == 4: 
                unmapped += 1
                continue
        
            try:
                if bwa_alg == 'bwasw' or use_stampy == 1:
                    #The SAM files output by the bwasw algorithm (or stampy) don't include the X0 or X1 tags
                    #bwasw created files do seem to include at least AS,XN,XS,XF,XE if we want to filter by those later.
                    ok = True
                else:
                    one_best_match = read.opt('X0') == 1
                    no_subpoptimal_matches = False
                    try: no_subpoptimal_matches = read.opt('X1') == 0
                    except: no_subpoptimal_matches = read.opt('XT')=='U'

                    ok = one_best_match and no_subpoptimal_matches

            except Exception, e: ## KeyError
                print '%s %s' % (read.qname, e)
                print read
                ok = False                    

            if ok:
                outfile.write(read)
            else:
                omit += 1
            
        print 'Removed %d reads out of %d' % (omit, i)
        print 'Omitted %d unmapped reads out of %d' % (unmapped, i)

if __name__ == '__main__':
    try:
        SamFilter().run()
    except Exception, e:
        print '%s' % e
        sys.exit(2)
