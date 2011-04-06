#!/usr/bin/env python
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

    def main(self):
        infile = pysam.Samfile(self.options.infile, 'r')
        outfile = pysam.Samfile(self.options.outfile, 'wh', template=infile)
        #outfile = pysam.Samfile(self.options.outfile, 'w', template=infile)
        print('Filtering reads in %s to %s' % (infile, outfile))
        i = 0
        omit = 0
        for read in infile.fetch():
            #if i > 0 and not i % 1e5: print i
            i += 1
  
            if not (read.flag in [0,4,16] and read.flag in [0,4,16]):
                print 'read %d flag %d not in {0,4,16}: omitting' % (i, read.flag)
                continue
            if read.flag == 4: continue
        
            try:
                one_best_match = read.opt('X0') == 1
                no_subpoptimal_matches = read.opt('X1') == 0
                indels = read.opt('XO') > 0 or read.opt('XG') > 0 

                ok = one_best_match and no_subpoptimal_matches ## and not indels

            except Exception, e: ## KeyError
                print '%s' % e
                ok = False

            if ok:
                outfile.write(read)
            else:
                omit += 1
            
        print 'Removed %d reads out of %d' % (omit, i)
        

if __name__ == '__main__':
    try:
        SamFilter().run()
    except Exception, e:
        print '%s' % e
        sys.exit(2)
