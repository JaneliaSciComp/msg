#!/usr/bin/env python

""" 
This is a tool to run a read simulator (specific for MSG reads) and generate a simulated reads
dataset along with a sam file showing how the reads should be mapped.

Runs mappers and measures runtime, number mapped, and mapping accuracy.

Notes:

1 - Generate a simulated data set of ~1,000,000 reads (from the attached mel X chrom is fine)

2 - map with default parameters for bwa, stampy and last.

3 - Calculate # mapped and mapping accuracy (compare read position from .fa file with position mapped in .sam file) for all three.

4 - determine speed of all three mappers

mel_X.fa | perl ~/msg/tools/simulate_SRshortreads_MSG_div.pl 95 1000000 200 5000

"""
from __future__ import division
import gzip
import shutil
import optparse
import os
import re
import sys
from time import time
import subprocess
from collections import defaultdict

import pysam
from Bio import SeqIO

__author__ = "Greg Pinero"
__author_email__ = "gregpinero@gmail.com"

############ USER SETTINGS ##################

# limits for distance of fragments in simulation
LOW_SIZE_LIM = 200
HIGH_SIZE_LIM = 5000
SIM_READS_FILE = 'reads.fa'
SIM_READS_SAM_FILE = 'reads.sam'
REF_FILE = 'dmel_X.fa'
STAMPY_PREMAP_W_BWA = True
STAMPY_SUBS_RATE = .03
SIMULATOR_PATH = '/home/pinerog/msg/tools/simulate_SRshortreads_MSG_div.pl'
# aln should be faster (and more suited for short reads).  bwasw may be more accurate. 
BWA_ALG = 'bwasw'
OVERWRITE_MAPPED_FILES = True

##############################################

def trace(fn):
    """A decorator to time your functions"""
    from time import time
    import sys
    def trace_func(*args, **kwargs):
        print fn.__name__ + '...',
        sys.stdout.flush()
        beg = time()
        ret = fn(*args, **kwargs)
        tot = time() - beg
        print '%.3f' % tot
        return ret
    return trace_func

def get_answers_by_qname(sim_reads_sam_file):
    """Get a dictionary of Direction	Start	CIGAR	MDtag by ReadID (qname)."""
    answers_by_qname = {}
    reads_file = open(sim_reads_sam_file)
    reads_file.next() #skip header line
    for line in reads_file:
        id, dir, start, cigar, mdtag = line.strip().split('\t')
        answers_by_qname[id] = (dir, start, cigar, mdtag)
    reads_file.close()
    return answers_by_qname

class Mapper(object):
    """
    This is a base object with methods that all mappers must implement.
    To create a new mapper, subclass this class and override at least the methods below
    that raise NotImplementedError's.    
    """

    #Have do_map set these at beginning and end of run
    def __init__(self, ref_file, reads_file):
        self.beg_time = self.end_time = None
        self.ref_file = ref_file
        self.reads_file = reads_file
        self.stdout_path = "%s.cmd.stdout" % self.__class__.__name__
        self.stderr_path = "%s.cmd.stderr" % self.__class__.__name__
        if os.path.exists(self.stdout_path): os.remove(self.stdout_path)
        if os.path.exists(self.stderr_path): os.remove(self.stderr_path)

    @trace
    def run_cmd(self, cmd):
        print cmd
        stdout = open(self.stdout_path,'a')
        stderr = open(self.stderr_path,'a')
        subprocess.check_call(cmd, shell=True, stdout=stdout, stderr=stderr)
        stdout.close()
        stderr.close()
        
    def do_map(self):
        raise NotImplementedError

    @trace
    def count_num_mapped(self):
        assert os.path.exists(self.out_filtered_sam_file_name)
        count = 0
        for line in open(self.out_filtered_sam_file_name):
            if not line.startswith('@'):
                count += 1
        return count
        
    @trace
    def calc_accuracy(self, answers_by_qname):
        """In the filtered sam file, compare the mapped location to the reference sam file of 
        known good locations"""
        assert os.path.exists(self.out_filtered_sam_file_name)
        num_correct = 0
        infile = pysam.Samfile(self.out_filtered_sam_file_name, 'r')
        for read in infile.fetch():
            assert read.qname in answers_by_qname, "Read %s found in test SAM file but not in simulated" % read.qname
            dir, start, cigar, mdtag = answers_by_qname[read.qname]
            #Interestingly, pysam seems to subtract 1 in read.pos compared to the actual sam file
            #The simulator currently reports 1 bp less than the sam files so everything works out.
            if int(dir) == int(read.flag) and int(start) == int(read.pos):
                num_correct += 1
        return num_correct            

    def get_run_time(self):
        assert self.beg_time and self.end_time,"do_map method must set self.beg_time and self.end_time"
        return '%.3f' % (self.end_time - self.beg_time)

    @property
    def out_sam_file_name(self):
        return "out.%s.sam" % self.__class__.__name__

    @property
    def out_filtered_sam_file_name(self):
        return "out.filtered.%s.sam" % self.__class__.__name__

    @trace
    def filter_sam_file(self):
        """Filter out bad reads, or multiply mapped reads"""
        
        flags_to_consider = set([0,4,16]) #no information, #segment unmapped, #SEQ being reverse complemented

        #Pass 1, get duplicate lines to remove
        infile = pysam.Samfile(self.out_sam_file_name, 'r')
        count_by_qname = defaultdict(int) #int returns 0
        for read in infile.fetch():
            if read.flag in set([0,16]):
                count_by_qname[read.qname] += 1
        infile.close()
        
        #Pass 2, figure out what to write
        infile = pysam.Samfile(self.out_sam_file_name, 'r')        
        outfile = pysam.Samfile(self.out_filtered_sam_file_name, 'wh', template=infile)
        #outfile = pysam.Samfile(self.options.outfile, 'w', template=infile)
        print('Filtering reads in %s to %s' % (self.out_sam_file_name, self.out_filtered_sam_file_name))
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
            #Remove multiple matches either by tag, or by seeing more than one
            if count_by_qname[read.qname] > 1:
                #print "removed read %s with %s mapped locations" % (read.qname, count_by_qname[read.qname])
                ok = False
            elif read.opt('X0'):
                one_best_match = read.opt('X0') == 1
                no_subpoptimal_matches = False
                try: no_subpoptimal_matches = read.opt('X1') == 0
                except: no_subpoptimal_matches = read.opt('XT')=='U'

                ok = one_best_match and no_subpoptimal_matches
            else:
                #if not multiple matches
                ok = True
            if ok:
                outfile.write(read)
            else:
                omit += 1
        print 'Removed %d reads out of %d' % (omit, i)
        print 'Omitted %d unmapped reads out of %d' % (unmapped, i)
        
        
class Stampy(Mapper):

    def pre_map_cleanup(self):
        subprocess.check_call('rm -f *.sthash *.stidx *.Stampy.sam', shell=True)

    def do_map(self):
        self.beg_time = time()
        self.pre_map_cleanup()
        self.run_cmd("bwa index -a bwtsw %s" % self.ref_file)
        self.run_cmd("stampy.py -G benchmark.stampy %s" % self.ref_file)
        self.run_cmd("stampy.py -g benchmark.stampy -H benchmark.stampy")
        if STAMPY_PREMAP_W_BWA:
            bwa_options = '--bwaoptions="-q10 %s"' % self.ref_file
        else:
            bwa_options = ''
        self.run_cmd("stampy.py %s -g benchmark.stampy -h benchmark.stampy --substitutionrate %s -M %s -o %s" % (
            bwa_options, STAMPY_SUBS_RATE, self.reads_file, self.out_sam_file_name))
        self.end_time = time()
    
class BWA(Mapper):
    def do_map(self):
        self.beg_time = time()
        self.run_cmd("bwa index -a bwtsw %s" % self.ref_file)        
        if BWA_ALG == 'aln':
            self.run_cmd("bwa %s -t 8 %s %s > bwa.sai" % (BWA_ALG, self.ref_file, self.reads_file))
            self.run_cmd("bwa samse %s bwa.sai %s > %s" % (self.ref_file, self.reads_file, self.out_sam_file_name))
        else:
            self.run_cmd("bwa %s -t 8 %s %s > %s" % (BWA_ALG, self.ref_file, self.reads_file,
            self.out_sam_file_name))     
        self.end_time = time()

class Last(Mapper):
    def do_map(self):
        """
        lastal options:
        -Q1 for fastq sanger
        -Q3 for fastq illumina

        lastdb options:
        -c ignore lowercase bases in ref?
        -m1111110 better for short strong alignments
        Can also set scoring:
        -r6 -q18 -a21 -b9 -e180 is match score, mismatch cost, gap cost of 21 + 9 * gap length, min aln score 180     
        is default of -e40 too high for short 95 bp reads?
        """
        self.beg_time = time()
        self.run_cmd("lastdb -c -v lastidx %s" % (self.ref_file))
        self.run_cmd("lastal -v lastidx %s > last.out.maf" % (self.reads_file))
        #Not sure if cmd below is usful??
        #self.run_cmd("cat last.out.maf | last-map-probs.py -s150 > last.out.filtered.maf")
        self.run_cmd("maf-convert.py sam last.out.maf > %s" % (self.out_sam_file_name))
        #Put in sq lines since last doesn't include it
        self.fix_headers()
        self.end_time = time()

    @trace
    def fix_headers(self):
        sq_lines = self.get_sq_header_lines()
        fixed_sam = open('last.sam.tmp','wb')
        fixed_sam.write(sq_lines + '\n')
        fixed_sam.write(open(self.out_sam_file_name).read())
        fixed_sam.close()
        shutil.move('last.sam.tmp', self.out_sam_file_name)    
    @trace
    def get_sq_header_lines(self):
        """Last doesn't include this due to a bug/oversight, so figure it out from ref file.
        Example: @SQ	SN:X	LN:22422827
        """
        #get order of ref id's in sam file
        sq_lines = []
        ids, seen = [], set()
        for line in open(self.out_sam_file_name):
            if not line.startswith('@'):
                ref_id = line.split('\t')[2].strip()
                if not ref_id in seen:
                    ids.append(ref_id)
                    seen.add(ref_id)
        assert ids
        #get ref file information and write out lines in the same order found in sam file above
        ref_io = SeqIO.to_dict(SeqIO.parse(open(self.ref_file),'fasta'))
        for ref_id in ids:
            sq_line = '\t'.join(["@SQ",'SN:'+ref_id,'LN:'+str(len(ref_io[ref_id].seq))])
            sq_lines.append(sq_line)
        return '\n'.join(sq_lines)

@trace
def run_simulator(options):
    if os.path.exists(SIM_READS_FILE) and os.path.exists(SIM_READS_SAM_FILE):
        print "Files already exist: %s, %s. Delete these to rerun simluator." % (SIM_READS_FILE,SIM_READS_SAM_FILE)
    else:
        cmd = "cat %s | perl %s %s %s %s %s" % (REF_FILE, SIMULATOR_PATH, options.maxlen, 
            options.numreads, LOW_SIZE_LIM, HIGH_SIZE_LIM)
        print "Running simulator:"
        print cmd
        subprocess.check_call(cmd, shell=True, stdout=open("sim.run.out",'w'),
            stderr=open("sim.run.err",'w'))
        assert os.path.exists(SIM_READS_FILE)
        assert os.path.exists(SIM_READS_SAM_FILE)

#### Set up Mappers to run
mappers = (Last, BWA, Stampy)

def main (argv=None):
    print "mapping benchmark starting run"
    
    if argv is None:
        argv = sys.argv
    
    usage = "Usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage, description=globals()['__doc__'])
    #parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
    
    parser.add_option ('--ref', help='Genome reference file to simulate reads from')
    parser.add_option ('--numreads', type='int', default=1000000, help='Number of reads to simulate')
    parser.add_option ('--maxlen', type='int', default=95, help='Maximum read length')
    try:
        (options, args) = parser.parse_args(argv[1:])
        if (not options.ref):
            parser.error('Must specify a genome reference file')
    except SystemExit: # Prevent exit when calling as function
        return 2

    #Run simulator
    run_simulator(options)
    answers = get_answers_by_qname(SIM_READS_SAM_FILE)

    results = ''

    #Run mappers and evaluate results
    for mapper in mappers:
        map = mapper(REF_FILE, SIM_READS_FILE)
        if os.path.exists(map.out_sam_file_name) and not OVERWRITE_MAPPED_FILES:
            run_time = 0
        else:
            map.do_map()
            map.filter_sam_file()
            run_time = map.get_run_time()
        num_mapped = map.count_num_mapped()
        num_correct = map.calc_accuracy(answers)
        results += """
Algorithm '%s' results:
  Mapping time: %s
  Number mapped: %s
  Number Correct: %s
""" % (map.__class__.__name__, run_time, num_mapped, num_correct)
    print results

if __name__ == '__main__':
    sys.exit(main())
    