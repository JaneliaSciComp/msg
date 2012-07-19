#!/usr/bin/env python

""" 
This is a tool to run a read simulator (specific for MSG reads) and generate a simulated reads
dataset along with a sam file showing how the reads should be mapped and then run the specified
mappers on the simulated file to compare runtime, number mapped, and mapping accuracy.


Command options:    
--ref (Genome reference file to simulate reads from)
--numreads (Number of reads to simulate)
--maxlen (Maximum read length)

Usage:
Example run:
(Generate simulated reads from the reference file dmel_X.fa with 1,000,000 reads, each
read having a maximum length of 200 and report the performance of the mappers)
python mapping_benchmark.py --ref dmel_X.fa --numreads 1000000 --maxlen 200

Also see the USER SETTINGS section of this file for more advanced options.

Adding a new mapper:

Create a new class (perhaps copy and existing one like Stampy for an easy start).  Make sure to
subclass the Mapper class.  Then simply code your own do_map method to have your mapper generate a sam
file.
Your do_map method should use these variables for the ref file, reads file: self.ref_file, self.reads_file
and output the final sam file to self.out_sam_file_name.
Make sure your do_map method has this statement as its first line: 
        self.beg_time = time()
and this statement as its last line:
        self.end_time = time()

Finally find the section in this code called "Set up Mappers to run" and add in your new mapper class.


Greg Pinero 2012 (gregpinero@gmail.com)
"""
from __future__ import division
import gzip
import shutil
import optparse
import os
import re
import sys
from time import time
import datetime
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
STAMPY_PREMAP_W_BWA = False
STAMPY_SUBS_RATE = .03
SIMULATOR_PATH = '/usr/local/msg/tools/simulate_SRshortreads_MSG_div.pl'
# aln should be faster (and more suited for short reads).  bwasw may be more accurate. 
BWA_ALG = 'bwasw'
OVERWRITE_MAPPED_FILES = True

##############################################

def get_free_memory():
    """
    Try to figure out how much memory is free on a Unix system.
    Returns free memory in mB.
    """ 
    data = open("/proc/meminfo", 'rt').readlines()
    free = 0
    for line in data:
        if line.startswith("MemFree") or line.startswith("Buffers") or line.startswith("Cached"):
            items = line.split()
            free += int(items[1])
    #print "free memory is",free/1000,"mB"
    return free/1000

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
        print '%.3f free: %s' % (tot, get_free_memory())
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
    
    See the adding a new mapper section at the top of this file for more information.
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
        """See the adding a new mapper section at the top of this file for more information
        on implementing this method."""
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
            if int(dir) == int(read.flag) and int(start) == int(read.pos)+1:
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
        self.run_cmd("maf-convert.py -d sam last.out.maf > %s" % (self.out_sam_file_name)) #-d to include @SQ header line
        self.end_time = time()

    def add_sq_header(self, file_name):
        """Last won't always include SQ headers, add them in here.
        WARNING: IF there is more than one chromosome I'm not sure they will be put in the right order, 
        fix it if it's a problem"""
        #@SQ	SN:X	LN:22422827
        headers=[]
        dref = SeqIO.to_dict(SeqIO.parse(open(self.ref_file),"fasta"))
        for key, val in dref.items():
            headers.append('\t'.join(['@SQ','SN:%s' % key, 'LN:%s' % len(val)]))
        new = open('last.tmp','w')
        for header in headers:
            new.write(header + '\n')
        new.write(open(file_name).read())
        new.close()
        shutil.move('last.tmp', file_name)

class LastForOptimize(Last):
    """ """
    def do_map(self, custom_args):
        self.beg_time = time()
        #self.run_cmd("lastdb -c -v lastidx %s" % (self.ref_file))
        self.run_cmd("lastal -v %s lastidx %s > last.out.maf" % (custom_args, self.reads_file))
        self.run_cmd("maf-convert.py sam last.out.maf > %s" % (self.out_sam_file_name)) #-d to include @SQ header line
        self.add_sq_header(self.out_sam_file_name)
        self.end_time = time()

class LastCustom2(Last):
    """ """
    def do_map(self):
        self.beg_time = time()
        #self.run_cmd("lastdb -c -v lastidx %s" % (self.ref_file))
        #self.run_cmd("lastal -v -r 259 -q 175 -a 24 -b 69 -x 99 -y 48 -z 180 -d 136 -e 285 -m 239 -l 45 lastidx %s > last.out.maf" % (self.reads_file))
        #Below got a really good score for 200bp
        self.run_cmd("lastal -v -r 263 -q 141 -a 237 -b 85 -x 68 -y 77 -z 112 -d 57 -e 205 -m 59 -l 35 lastidx %s > last.out.maf" % (self.reads_file))
        self.run_cmd("maf-convert.py sam last.out.maf > %s" % (self.out_sam_file_name)) #-d to include @SQ header line
        self.add_sq_header(self.out_sam_file_name)
        self.end_time = time()

class LastCustom(Last):
    """ """
    def do_map(self):
        self.beg_time = time()
        self.run_cmd("lastdb -c -v lastidx %s" % (self.ref_file))
        self.run_cmd("lastal -v lastidx %s | last-map-probs.py -s20 > lastbd.out.maf" % (self.reads_file))
        #self.run_cmd("lastal -v -r 263 -q 141 -a 237 -b 85 -x 68 -y 77 -z 112 -d 57 -e 205 -m 59 -l 35 lastidx %s | last-map-probs.py -s20 > lastbd.out.maf" % (self.reads_file))
        self.run_cmd("maf-convert.py sam lastbd.out.maf > %s" % (self.out_sam_file_name)) 
        #-d didn't help here to include @SQ header line
        self.add_sq_header(self.out_sam_file_name)
        self.end_time = time()

class LastBestProbs(Last):
    """Test Last with filtering to keep best prob. alignments"""
    def do_map(self):
        self.beg_time = time()
        self.run_cmd("lastdb -c -v lastidx %s" % (self.ref_file))
        self.run_cmd("lastal -v lastidx %s | last-map-probs.py -s150 > lastb.out.maf" % (self.reads_file))
        self.run_cmd("maf-convert.py sam lastb.out.maf > %s" % (self.out_sam_file_name)) 
        #-d didn't help here to include @SQ header line
        self.add_sq_header(self.out_sam_file_name)
        self.end_time = time()

@trace
def run_simulator(options):
    if os.path.exists(SIM_READS_FILE) and os.path.exists(SIM_READS_SAM_FILE):
        print "Files already exist: %s, %s. Delete these to rerun simluator." % (SIM_READS_FILE,SIM_READS_SAM_FILE)
    else:
        cmd = "cat %s | perl %s %s %s %s %s" % (options.ref, SIMULATOR_PATH, options.maxlen, 
            options.numreads, LOW_SIZE_LIM, HIGH_SIZE_LIM)
        print "Running simulator:"
        print cmd
        subprocess.check_call(cmd, shell=True, stdout=open("sim.run.out",'w'),
            stderr=open("sim.run.err",'w'))
        assert os.path.exists(SIM_READS_FILE)
        assert os.path.exists(SIM_READS_SAM_FILE)

#### Set up Mappers to run #########
mappers = (LastBestProbs, LastCustom, Last)

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
        map = mapper(options.ref, SIM_READS_FILE)
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
    open(datetime.datetime.today().isoformat().replace(':',"_")+'.results','w').write(results)

if __name__ == '__main__':
    sys.exit(main())
    