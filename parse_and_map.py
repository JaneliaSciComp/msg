#!/usr/bin/env python

##parse and map
##takes Illumina file and barcodes file and parses into individual files
##then maps reads using bwa to towo parental genomes
##also exports stats on parsed files
##0.2.1 provide option to just map pre-parsed reads
##0.2.2 Now uses correct Andolfatto parsing code (parse_BCdata2BWA.v9_laptop.pl). Fixes error that was discarding many good reads as bad_barcodes.
##      Provide options for which genome sequence to map to.
##0.2.3 Allow user to just parse data
##0.2.4 multithread for 4 processors, added option for sec_velvet_61 genome
##0.2.5 have used provide full name of genome
##0.2.6 add option to map a subset of individuals
##0.3.0 Dan: use CommandLineApp class with ability to take all input parameters on command line
__version__ = '0.3.0'

import subprocess
import sys, os
import glob
import time, datetime
from mapping_functions import *
from msglib import *
from cmdline.cmdline import CommandLineApp

##### INTERNAL OPTIONS (for developers) ######

# Should MD tags be added to mapped SAM files when they are not included by default.
# It seems to be very slow to add these in, but it does affect the end result.  I'm not sure 
# where though.
GEN_MD = True

# Be verbose
DEBUG = True

#############################

class ParseAndMap(CommandLineApp):
    def __init__(self):
        CommandLineApp.__init__(self)
        
        op = self.option_parser
        op.set_usage('usage: parse_and_map -[npm] -i reads.fq -b barcodes --parent1 parent1_genome --parent2 parent2_genome')

        op.add_option('-i', '--raw-data', dest='raw_data_file', type='string', default=None, 
                      help='Raw Illumina data file')

        op.add_option('-b', '--barcodes', dest='barcodes_file', type='string', default=None, 
                      help='Barcodes file')

        op.add_option('--re_cutter', dest='re_cutter', type='string', default=None, 
                      help='restriction enzyme')

        op.add_option('--linker_system', dest='linker_system', type='string', default=None, 
                      help='linker system')

        op.add_option('--parent1', dest='parent1', type='string', default=None, 
                      help='Parent 1 genome file to align against')

        op.add_option('--parent2', dest='parent2', type='string', default=None, 
                      help='Parent 2 genome file to align against')

        op.add_option('-m', '--map-only', dest='map_only', default=False, action='store_true', 
                      help='Just map data?')

        op.add_option('-n', '--num-ind', dest='num_ind', default=None, type='int',
                      help='Number of individuals to map')

        op.add_option('-p', '--parse-only', dest='parse_only', default=False, action='store_true',
                      help='Just parse data?')

        op.add_option('--bwa_alg', dest='bwa_alg', type='string', default="aln", 
                      help='Algorithm for BWA mapping. Use aln or bwasw. This is ignored is use_stampy is 1.')

        op.add_option('--bwa_threads', dest='bwa_threads', type='int', default=1, 
                      help='Number of threads in the BWA multi-threading mode')

        op.add_option('--use_stampy', dest='use_stampy', type='int', default=0, 
                      help='Use stampy for mapping')

        op.add_option('--stampy_premap_w_bwa', dest='stampy_premap_w_bwa', type='int', default=1, 
                      help='Use BWA before running stampy for a speed boost.')

        op.add_option('--stampy_pseudo_threads', dest='stampy_pseudo_threads', type='int', default=0, 
                      help='Use qsub commands to split up stampy processes during mapping. NOT currently used.')

        #Set divergence for mapping to a foreign reference (note: this is
        #strongly recommended for divergences >3%, as it will automatically
        #shut down BWA pre-mapping which otherwise cause occasional segfaults)
        #(stampy default if .001)
        op.add_option('--indiv_stampy_substitution_rate', dest='stampy_substitution_rate', type='float', default=0.001, 
                      help='Set divergence for mapping to a foreign reference')
                      
        op.add_option('--indiv_mapq_filter', dest='mapq_filter', type='int', default=0, 
                      help='Filter out poor alignments.  Set this to 0 to skip.')
                      

    def main(self):
        print datetimenow()
        print bcolors.OKBLUE + 'parse_and_map version %s' % __version__ + bcolors.ENDC

        self.start_time = time.time()

        # Get variables from barcode file
        #self.variables, self.bc = read_barcodes(self.options.barcodes_file)
        self.bc = read_barcodes(self.options.barcodes_file)
        
        # determine number to map 0.2.6
        if not self.options.num_ind: # if no entry, default is all individuals
            self.options.num_ind = len(self.bc)
        elif int(self.options.num_ind) > len(self.bc): # if user accidentally enters more than the number of individuals, run # individuals
            self.options.num_ind = len(self.bc)
        elif self.options.num_ind.lower() == "all":
            self.options.num_ind = len(self.bc) # if user enters "[Aa]ll"

        self.logdir = 'log'
        if not os.path.exists(self.logdir):
            os.mkdir(self.logdir)

        self.parsedir = os.path.basename(self.options.raw_data_file) + '_parsed'
        self.samdir = os.path.basename(self.options.raw_data_file) + '_sam_files'

        self.parsed_time = time.time()
        if self.options.map_only:
            print bcolors.WARN + 'Refusing to parse raw reads: --map-only option is in effect' + bcolors.ENDC

        elif os.path.exists(self.parsedir):
            print bcolors.WARN + 'Refusing to parse raw reads: parsed output directory %s exists' % self.parsedir + bcolors.ENDC
        else:
            self.parse()

        if self.options.parse_only:
            print bcolors.WARN + 'Refusing to map parsed reads: --parse-only option is in effect' + bcolors.ENDC
#        elif os.path.exists(self.samdir):
#            print bcolors.WARN + 'Refusing to map parsed reads: samfiles output directory %s exists' % self.samdir + bcolors.ENDC
        else:
            self.map()
            if self.options.bwa_alg == 'aln':
                #bwasw and stampy don't make sai files    
                self.delete_files()

    def parse(self):

        # Convert Illumina reads to format for BWA - Peters program does not create a new file,
        # instead it writes to existing files. So first need to remove the older files.
        # check if output files already exist. If so, erase them.
        if os.path.isfile("./bad_barcodes"): # shorthand-if bad_barcodes file exists, then so do others
            os.remove("./bad_barcodes")
            os.remove("./linker")
            for ind in self.bc:
                fastq_file = 'indiv' + ind[1] + '_' + ind[0]
                os.remove("./" + fastq_file)
    
        print "Parsing data into individual barcode files"

    ## USAGE: parse_BCdata2BWA.pl version 11
    ##    fastq_file
    ##    dir_name
    ##    enzyme(MseI/NdeI/Hpy188I/0=null)
    ##    prepend(0|1)
    ##    barcode_file
    ##    ignored_barcode_file(null)
    ##    linker_type(Dros_SR_vI|Dros_SR_vII|Dros_SR_vIII_T_overhang|Dros_SR_vIII_blunt|ButtFinch_PE_vI|ButtFinch_PE_vII)
    ##    convert_QV(0|1)
    ##    strict(0|1)
    ##    minQV(default=0)
    ##    assign barcode if first n-1 bps are nonredundant(default=0)
    ##     when strict = 0,
    ##         junk = sequences < 20bp and sequences with mononucleotide runs >8
    ##      when strict = 1,
    ##          the sequences are also checked to all start with the correct bases (corresponding to the RE used)
    ##          and dumped to junk if they don't start with the right base(s)
    ##
    ##  The following should be provided in the barcodes file
    ##  sp1  sp2  enzyme  prepend  linker_type  convert_QV  strict  minQV

    ## OPTIONS: 04.05.11
    ##    -e restriction enzyme used for digestion [MseI|NdeI|Hpy188I|null]
    ##    -p prepend bases to reads to complete the RE motif (useful for mapping reads) [0|1]
    ##    -b barcode file [barcode_file]
    ##    -i file of barcodes to ignore [null]
    ##    -l linker system used [Dros_SR_vI|Dros_SR_vII|Dros_SR_vIII_T_overhang|Dros_SR_vIII_blunt|ButtFinch_PE_vI|ButtFinch_PE_vII]
    ##    -s strict (only retain reads with the correct RE motif) [0|1]
    ##    -m minimum QV to mask from the 3' end [0]
    ##    -a assign barcode if first n-1 bps are nonredundant [0]
    ##    -c convert individual fastq files into fasta [0]
    
        raw_data = self.options.raw_data_file
        args = ["perl", os.path.join(os.path.dirname(__file__), "parse_BCdata2BWA.pl"), 
                  '-b', self.options.barcodes_file,
                  '-e', self.options.re_cutter,
                  '-l', self.options.linker_system,
                  raw_data ]

        print ' '.join(args)
        sys.stdout.flush()
        subprocess.call(args) 


        # Outputs separate file for each individual with barcode and identifier in filename
        # "indiv#_barcode"
    
    
        #print file with statistics of parsing
        stats_file = open(raw_data+'_stats.txt','w')
        total_reads = 0
        read_file = './'+raw_data+'_parsed/bad_barcodes'
        number_reads = count_lines(read_file)/4
        stats_file.write('bad_barcodes\t%s\n' %(number_reads))
        total_reads += number_reads
    
        read_file = './'+raw_data+'_parsed/unreadable_barcodes'
        number_reads = count_lines(read_file)/4
        stats_file.write('unreadable_barcodes\t%s\n' %(number_reads))   
        total_reads += number_reads
    
        read_file = './'+raw_data+'_parsed/linkers'
        number_reads = count_lines(read_file)/2
        stats_file.write('linkers\t%s\n' %(number_reads))   
        total_reads += number_reads
    
        read_file = './'+raw_data+'_parsed/junk'
        number_reads = count_lines(read_file)/4
        stats_file.write('junk\t%s\n' %(number_reads))   
        total_reads += number_reads

        for ind in self.bc:
            file_name = 'indiv' + ind[1] + '_' + ind[0]
            fastq_file = raw_data + '_parsed/' + file_name
            number_reads = count_lines(fastq_file)/4
            total_reads += number_reads
            stats_file.write('indiv%s_%s\t%s\n' %(ind[1],ind[0],number_reads))
            #TEMP - make sym links to parsed fq files with .fq extensions to they can run in stampy
            #Later versions of stampy will hopefully fix this so we can remove it.
            #ln trivia: symlink target should be relative to sym link, not where you are.
            sym_link_args = 'ln -s %s %s.fq' % (file_name, fastq_file)
            print "create ln with .fq extension:", sym_link_args
            subprocess.check_call([sym_link_args], shell=True)
                
        stats_file.write('total_reads\t%s' %(total_reads))
        stats_file.close()
    
        self.parsed_time = time.time()
        parsing_time = (self.parsed_time - self.start_time)
        print "Parsing took about %s minutes" %(parsing_time/60)
    
    def _map_w_stampy(self, fastq_file, parent1, parent2, aln_par1_sam, aln_par2_sam, file_par1_log,
        file_par2_log, misc_indiv_log):
        #Align each to sim - output fastq file
        #TEMP: stampy requires .fq extension on our files.  Remove this if stampy fixes it.
        # I'm temporarily creating symlinks in the parse step.
        fastq_file += '.fq'
        if self.options.stampy_premap_w_bwa == 1:
            args = ['stampy.py', '-v3', '--inputformat=fastq', '--substitutionrate=%s' % self.options.stampy_substitution_rate,
                '--bwaoptions="-q10 %s"' % parent1, '-g',  
                "%s.stampy.msg" % parent1, '-h', "%s.stampy.msg" % parent1, "-M",
                fastq_file, "-o", aln_par1_sam]
        else:
            args = ['stampy.py', '--inputformat=fastq', '-g', "%s.stampy.msg" % parent1, 
                '--substitutionrate=%s' % self.options.stampy_substitution_rate,
                '-h', "%s.stampy.msg" % parent1, "-M",
                fastq_file, "-o", aln_par1_sam]
        print "stampy call parent 1:"
        print ' '.join(args)
        #popen notes: When shell==True, send in only one argument in list
        file_par1_sam = subprocess.Popen([' '.join(args)],
            stderr=file_par1_log, shell=True)
        
        #Align each to sec - output fastq file
        if self.options.stampy_premap_w_bwa == 1:
            args = ['stampy.py', '-v3', '--inputformat=fastq', '--substitutionrate=%s' % self.options.stampy_substitution_rate,
                '--bwaoptions="-q10 %s"' % parent2, '-g',  
                "%s.stampy.msg" % parent2, '-h', "%s.stampy.msg" % parent2, "-M",
                fastq_file, "-o", aln_par2_sam]
        else:
            args = ['stampy.py', '--inputformat=fastq', '-g', "%s.stampy.msg" % parent2, 
                '--substitutionrate=%s' % self.options.stampy_substitution_rate,
                '-h', "%s.stampy.msg" % parent2, "-M",
                fastq_file, "-o", aln_par2_sam]
        print "stampy call parent 2:"
        print ' '.join(args)
        file_par2_sam = subprocess.Popen([' '.join(args)],
            stderr=file_par2_log, shell=True)

        #pause until these two processes are finished. This is a precaution. Don't want to continue until sure sam files are completely written 
        file_par1_sam.wait()
        file_par2_sam.wait()
        
        #Fix stampy generated sam files.  Fix headers, and sort
        for file_to_fix in (aln_par1_sam, aln_par2_sam):
            #Remove @PG header line from STAMPY generated SAM files since PYSAM dies on the PN: field.
            subprocess.check_call(['grep -v "@PG" %s > %s.fixed' % (file_to_fix, file_to_fix)], 
                shell=True, stdout=misc_indiv_log, stderr=misc_indiv_log)
            subprocess.check_call(['mv %s.fixed %s' % (file_to_fix, file_to_fix)], shell=True,
                stdout=misc_indiv_log, stderr=misc_indiv_log)    
            #convert to bam to prepare for sorting  
            subprocess.check_call(['samtools view -btSh -o %s.bam %s' % (file_to_fix, file_to_fix)],
                shell=True, stdout=misc_indiv_log, stderr=misc_indiv_log)
            #samtools adds .bam suffix to output FYI
            subprocess.check_call(['samtools sort -n %s.bam %s.sorted' % (file_to_fix, file_to_fix)],
                shell=True, stdout=misc_indiv_log, stderr=misc_indiv_log)
            #convert back to SAM
            subprocess.check_call(['samtools view -h -o %s %s.sorted.bam' % (file_to_fix, file_to_fix)],
                shell=True, stdout=misc_indiv_log, stderr=misc_indiv_log)
            #remove temporary sorting files
            os.remove('%s.bam' % file_to_fix)
            os.remove('%s.sorted.bam' % file_to_fix)
            
    def map(self):

        #Run bwa programs on each indiv file
        #Output to new folder
    
        raw_data = os.path.basename(self.options.raw_data_file)
        #Make directory to hold sam files
        dirname = self.samdir
        if not os.path.isdir("./" + dirname):
            os.mkdir("./" + dirname)

        par1 = 'par1'
        par2 = 'par2'
        parent1 = self.options.parent1
        parent2 = self.options.parent2

        print bcolors.OKBLUE + "Mapping reads to genomes" + bcolors.ENDC
    
        barcodes_file = open(self.options.barcodes_file,'r')
        barcodes_file.readline()##ignore first two lines of barcodes file
        barcodes_file.readline()
        sample_num = 0 ##0.2.6
        for ind in self.bc:
            sample_num +=1 ##0.2.6
            fastq_file_name = 'indiv' + ind[1] + '_' + ind[0]
            fastq_file = './' + raw_data + '_parsed/' + fastq_file_name

            #Change format for sim - output sam file
            aln_par1_sam = './' + raw_data + '_sam_files/aln_' + fastq_file_name + "_" + par1 + ".sam"
            #Change format for sec - output sam file
            aln_par2_sam = './' + raw_data + '_sam_files/aln_' + fastq_file_name + "_" + par2 + ".sam"

            if ((os.path.exists(aln_par1_sam) or os.path.exists(aln_par1_sam+'.gz')) and
                (os.path.exists(aln_par2_sam) or os.path.exists(aln_par2_sam+'.gz'))): 
                print bcolors.WARN + 'Refusing to map reads for ' + fastq_file_name + bcolors.ENDC
                continue

            file_par1_log = open(os.path.join(self.logdir, fastq_file_name + par1 + '.log'), "w")
            file_par2_log = open(os.path.join(self.logdir, fastq_file_name + par2 + '.log'), "w")
            misc_indiv_log = open(os.path.join(self.logdir, fastq_file_name + '.misc.log'), "w")

            if self.options.use_stampy == 1:
                self._map_w_stampy(fastq_file, parent1, parent2, aln_par1_sam, aln_par2_sam, file_par1_log,
                    file_par2_log, misc_indiv_log)

            elif self.options.bwa_alg == 'aln':
                #Align each to sim - output fastq file
                aln_par1_sai =  './aln_' + fastq_file_name + "_" + par1 + ".sai"
                file_par1_sai = open(aln_par1_sai,"w")
                file_par1_sai = subprocess.Popen(['bwa', 'aln', 
                    '-t ' + str(self.options.bwa_threads), parent1, fastq_file],
                    stdout=file_par1_sai, stderr=file_par1_log)
                
                #Align each to sec - output fastq file
                aln_par2_sai =  './aln_' + fastq_file_name + "_" + par2 + ".sai"
                file_par2_sai = open(aln_par2_sai,"w")
                file_par2_sai = subprocess.Popen(['bwa', 'aln', '-t ' + str(self.options.bwa_threads), 
                    parent2, fastq_file], 
                    stdout=file_par2_sai, stderr=file_par2_log)

                #pause until these two processes are finished. If don't, then samse starts on empty or incomplete file
                file_par1_sai.wait()
                file_par2_sai.wait()

                file_par1_sam = open(aln_par1_sam,'w')
                file_par1_sam = subprocess.Popen(['bwa', "samse", parent1, aln_par1_sai, fastq_file],stdout=file_par1_sam)
                file_par2_sam = open(aln_par2_sam,'w')
                file_par2_sam = subprocess.Popen(['bwa', "samse", parent2, aln_par2_sai, fastq_file],stdout=file_par2_sam)

                #pause until these two processes are finished. This is a precaution. Don't want to continue until sure sam files are completely written 
                file_par1_sam.wait()
                file_par2_sam.wait()

            elif self.options.bwa_alg == 'bwasw':
                #Align each to sim - output fastq file
                file_par1_sam = open(aln_par1_sam,'w')
                file_par1_sam = subprocess.Popen(['bwa', 
                    'bwasw', '-t ' + str(self.options.bwa_threads), parent1, fastq_file],
                    stdout=file_par1_sam, stderr=file_par1_log)
                
                #Align each to sec - output fastq file
                file_par2_sam = open(aln_par2_sam,'w')
                file_par2_sam = subprocess.Popen(['bwa', 'bwasw', 
                    '-t ' + str(self.options.bwa_threads), parent2, fastq_file], 
                    stdout=file_par2_sam, stderr=file_par2_log)

                #pause until these two processes are finished. This is a precaution. Don't want to continue until sure sam files are completely written 
                file_par1_sam.wait()
                file_par2_sam.wait()

            else:
                raise ValueError("Not using stampy and invalid bwa_alg option: %s. Use aln or bwasw" % self.options.bwa_alg)
   
            #After updating files with options below, should we keep the intermediate version around:
            put_back_command = DEBUG and 'cp' or 'mv' #means 'cp' if DEBUG else 'mv'
            if self.options.mapq_filter:
                # remove poor alignments if requested
                for (sam_file, log_file) in ((aln_par1_sam,file_par1_log), (aln_par2_sam,file_par2_log)):
                    subprocess.check_call('samtools view -Sh -q %s -o %s.mapq_filtered.sam %s' % (
                        self.options.mapq_filter ,sam_file, sam_file),
                        shell=True, stdout=log_file, stderr=log_file)
                    result = subprocess.check_call("%s -f %s %s" % (put_back_command, sam_file + '.mapq_filtered.sam',sam_file), shell=True)
   
            if GEN_MD and (self.options.bwa_alg == "bwasw" or self.options.use_stampy == 1):
                #Add in MD tags since bwasw omits these
                #(Write out to <output>.tmp.sam and then move.  Don't overwite input file directly since piped commands outputs continually.
                #TODO: It might be worth sorting the input files first to speed this up. Measure and test.
                for (sam_file, parent_, log_file) in ((aln_par1_sam, parent1, file_par1_log),(aln_par2_sam, parent2, file_par2_log)):
                    result = subprocess.check_call(
                        "samtools calmd -uS %s %s | samtools view -h -o %s -" % (sam_file, parent_, sam_file + '.added_calmd.sam'), 
                        shell=True, stderr=log_file)
                    result = subprocess.check_call("%s -f %s %s" % (put_back_command, sam_file + '.added_calmd.sam',sam_file), shell=True)

            print "done sample " + fastq_file   

            if int(self.options.num_ind) == sample_num:##0.2.6
                break
            
        barcodes_file.close()    ##0.2.5
    
        self.mapped_time = time.time()
        mapping_time = self.mapped_time - self.parsed_time
        print "Mapping took about %s minutes" %(mapping_time/60)
            
    def delete_files(self):
        #Delete sai files
        print "Deleting .sai files"   ##rewritten for 0.2.5 to delete only files processed in this script

        barcodes_file = open(self.options.barcodes_file,'r')
        barcodes_file.readline()##ignore first two lines of barcodes file
        barcodes_file.readline()
        sample_num = 0##0.2.6
        for ind in self.bc:
            sample_num +=1##0.2.6
            fastq_file_name = 'indiv' + ind[1] + '_' + ind[0]
            target1 = './aln_' + fastq_file_name + '_par1.sai'
            target2 = './aln_' + fastq_file_name + '_par2.sai'

            if os.path.exists(target1): os.remove(target1)
            else: print bcolors.WARN + 'Missing .sai file: ' + target1 + bcolors.ENDC
            if os.path.exists(target2): os.remove(target2)
            else: print bcolors.WARN + 'Missing .sai file: ' + target2 + bcolors.ENDC

            if int(self.options.num_ind) == sample_num:##0.2.6
                break
        barcodes_file.close()  
    
    def genotype(self):
        #Make genomewide genotype calls for each individual
        #Output to new folder as EXCEL file
        print "Parsing bwa output into genomewide genotype calls"
    
    
        #Make directory to hold genotype calls
        dirname = self.options.raw_data_file + "_genotypes"
        if not os.path.isdir("./" + dirname + "/"):
            os.mkdir("./" + dirname + "/")
    
        barcodes_file = open(self.options.barcodes_file,'r')
        sample_num = 0##0.2.6
        for ind in self.bc:
            sample_num +=1##0.2.6
            fastq_file = 'indiv' + ind[1] + '_' + ind[0]
        ##        file1 = './' + raw_data + '_sam_files/aln_' + fastq_file + "_" + sp1 + ".sam"
        ##        file2 = './' + raw_data + '_sam_files/aln_' + fastq_file + "_" + sp2 + ".sam"
    
        #Peter change to output CORRECT file names to new directory ./genotypes
    
        #Peter change to raw_data genotype directoryd
    
            subprocess.call(["perl", "parse_BWA2sp.v8.3.pl", self.variables[0], self.variables[1], fastq_file, self.options.raw_data_file]) 
    
            if int(self.options.num_ind) == sample_num:##0.2.6
                break
    
        barcodes_file.close()
        genotyped_time = time.time()
        genotyping_time = genotyped_time - mapped_time
        print "Genotyping every marker took about %s minutes" %(genotyping_time/60)

    def ask_user_if_options_not_specified(self):
        if not self.options.raw_data_file:
            self.options.raw_data_file = raw_input("Enter the Illumina data filename: ")
            
        if not self.options.barcodes_file:
            self.options.barcodes_file = raw_input("Enter name of barcodes file (default = 'barcodes') ")
        if not self.options.barcodes_file:
            self.options.barcodes_file = 'barcodes'

        if not self.options.map_only and not self.options.parse_only:
            parse_or_map = raw_input("Parse and map data (default), just parse (P), or just map parsed reads (M)? ")##0.2.3
            parse_or_map = parse_or_map.lower()
            assert parse_or_map in ['','m','p']
            if parse_or_map == 'm':
                self.options.map_only = True
            elif parse_or_map.lower() == 'p':
                self.options.parse_only = True

        if not self.options.numind:
            self.options.num_ind = raw_input("How many individuals to map (default = All)? ")##0.2.6



def datetimenow():
    return str(datetime.datetime.now()).split('.')[0].replace(' ', '_').replace(':','.')

if __name__ == '__main__':
    try:
        ParseAndMap().run()
    except Exception, e:
        print 'Error in parse_and_map:\n'
        print '%s' % e
        sys.exit(2)
