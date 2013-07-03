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
import shutil
import glob
import time, datetime
import gzip
from mapping_functions import *
from msglib import *
from cmdline.cmdline import CommandLineApp
import barcode_splitter

##### INTERNAL OPTIONS (for developers) ######

# Should MD tags be added to mapped SAM files when they are not included by default.
# It seems to be very slow to add these in, but it does affect the end result.  I'm not sure 
# where though.
GEN_MD = True
GZIP_OUTPUT = True

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
        
        op.add_option('--dbg', dest='debug', default=False, action='store_true', 
                      help='More verbose output and leaves temporary files instead of cleaning up')

        op.add_option('--quality_trim_reads_thresh', dest='quality_trim_reads_thresh', default=None, type='int',
                      help='Illumina parsing only - what level should split files be quality trimmed at')
                      
        op.add_option('--quality_trim_reads_consec', dest='quality_trim_reads_consec', default=None, type='int',
                      help='Illumina parsing only - Minimum number of consecutive bases passing threshold values')

        #Set divergence for mapping to a foreign reference (note: this is
        #strongly recommended for divergences >3%, as it will automatically
        #shut down BWA pre-mapping which otherwise cause occasional segfaults)
        #(stampy default if .001)
        op.add_option('--indiv_stampy_substitution_rate', dest='stampy_substitution_rate', type='float', default=0.001, 
                      help='Set divergence for mapping to a foreign reference')
                      
        op.add_option('--indiv_mapq_filter', dest='mapq_filter', type='int', default=0, 
                      help='Filter out poor alignments.  Set this to 0 to skip.')

        op.add_option('--new_parser', dest='new_parser', type='int', default=0, 
                      help='Use a faster, experimental parser')
                
        #Illumina indexing only
        op.add_option('--index_file', dest='index_file', type='string', default=None, 
                      help='When using Illumina indexes, this is the fastq/fasta file with the indexes.')
        op.add_option('--index_barcodes', dest='index_barcodes', type='string', default=None, 
                      help='When using Illumina indexes, this the the file with a listing of index sequences and labels')
                      
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
            if self.options.index_file and self.options.index_barcodes:
                final_paths = self.parse_illumina_indexes()
            else:
                final_paths = self.parse()
            self.qual_trim_parsed_files(final_paths)

        if self.options.parse_only:
            print bcolors.WARN + 'Refusing to map parsed reads: --parse-only option is in effect' + bcolors.ENDC
#        elif os.path.exists(self.samdir):
#            print bcolors.WARN + 'Refusing to map parsed reads: samfiles output directory %s exists' % self.samdir + bcolors.ENDC
        else:
            self.map()
            if self.options.bwa_alg == 'aln':
                #bwasw and stampy don't make sai files    
                self.delete_files()

    def qual_trim_parsed_files(self, final_paths):
        """Optionally quality trim the files.
        Go through all parsed files and quality trim them.
        The process should be the same for both illumina indexes, and normal parsing.
        Maintain existing names"""
        if self.options.quality_trim_reads_thresh and self.options.quality_trim_reads_consec:
            for path in final_paths:
                if GZIP_OUTPUT:
                    gzip_switch = '-z'
                    new_expected_file_name = path + '.trim.fastq.gz'
                else:
                    gzip_switch = ''
                    new_expected_file_name = path + '.trim.fastq'
                args = [os.path.join(os.path.dirname(__file__), "TQSfastq.py"), 
                          '-f', path, '-t', str(self.options.quality_trim_reads_thresh),
                          '-c', str(self.options.quality_trim_reads_consec), '-q', gzip_switch,
                          '-o', path]
                print ' '.join(args)
                subprocess.check_call(' '.join(args), shell=True)
                if self.options.debug:
                    shutil.copy(path, path + '.pretrim')
                os.remove(path)
                shutil.move(new_expected_file_name, path)
            
    def parse_illumina_indexes(self):
        """
        When the illumina index system is used.  We have 3 input files: A normal fastq file with all of the reads,
        a fastq with the same reference ids as the reads file but with the index as the sequence, and a small
        indexes file that contains a label for each index sequence.
        
        This method calls barcodes_splitter.py which will parse out each index into the working directory.
        The strategy here is to grab the relevant split files output from barcodes_splitter.py and call the
        standard issue MSG parser on each of them
        and then move and rename those results to where/what MSG would expect.   Also note msgCluster.pl calls the 
        make_msg_barcodes_file function in barcode_splitter to create an updated barcodes file 
        at startup so downstream MSG programs know where to find the parsed data.
        """
        print "Starting illumina index parse"
        prefix = 'index_parse' #used to keep track of index parsing output files

        #store and return the final relative paths of all the parsed individuals' fastq files
        final_paths = []

        #capture output here
        log = open(os.path.join(self.logdir, '%s.log' % (prefix)), "w")

        def run_command(args):
            """simple boilerplate to run commands"""
            print ' '.join(args)
            #suproccess.check_call failed silently when this was run from qsub using a 
            #list as args, so I join it into a string which seems to work everywhere.
            subprocess.check_call(' '.join(args), shell=True, stdout=log, stderr=log)

        #Mkdir output directory
        if not os.path.exists(self.parsedir):
            os.mkdir(self.parsedir)

        #call barcode_splitter.py to break up by illumina indexes
        args = ['python', os.path.join(os.path.dirname(__file__), "barcode_splitter.py"), 
                  '--bcfile', str(self.options.index_barcodes),
                  '--prefix', prefix,
                  '--idxread 1', self.options.index_file, self.options.raw_data_file ]
        if self.options.index_file.endswith('.gz'):
            assert self.options.raw_data_file.endswith('.gz'), "Both of these files must be compressed or un-compressed: %s, %s" % (
                self.options.index_file, self.options.raw_data_file)
            #make sure output gets a .gz suffix
            args.insert(7, '--suffix')
            args.insert(8, '.gz')
        run_command(args)
        
        #Gather and process output files
        barcodes_dict = barcode_splitter.read_barcodes(self.options.index_barcodes) #id by seq
        for index_id in set(barcodes_dict.values()):
            #Ignore all of the *_1 files because they came from illumina index fastq file 
            expected_file_name = "%s%s_read_2" % (prefix, index_id)
            if self.options.index_file.endswith('.gz'):
                expected_file_name += '.gz'
            print "processing file",expected_file_name
            assert os.path.exists(expected_file_name)
            
            #Process file through regular MSG parser
            self.parse(expected_file_name)
            expected_msg_parse_dir = expected_file_name + '_parsed' #this is where it will put the parsed files
            
            #Copy parsed files to output directory and rename files to denote which index they came from
            def make_new_name(old_name, index_id):
                """example rename for an illumina index called "standard"
                indivA12_AATAAG -> indivA12standard_AATAAG
                Note that barcode_splitter.make_msg_barcodes_file expects the file 
                names to be in this format, so update that too if changing naming here.
                """
                if old_name.startswith('indiv'):
                    parts = old_name.split('_')
                    return parts[0] + index_id + '_' + parts[1]
                else:
                    return old_name + '_' + index_id
                    
            for fn in os.listdir(expected_msg_parse_dir):
                #Skip .fq files since we assume they are symlinks. (this gets hairy).  
                #Explanation: The MSG parser created symbolic links to each file appending an .fq so
                #stampy can process them.  Since we're renaming the files we need to not copy the 
                #symlinks and later recreate them
                if not fn.endswith('.fq'):
                    #example rename indivA12_AATAAG -> indivA12standard_AATAAG
                    new_name = make_new_name(fn, index_id)
                    shutil.copy(os.path.join(expected_msg_parse_dir, fn), os.path.join(self.parsedir, new_name))
                    #Recreate symlink we didn't copy and store final path
                    if fn.startswith('indiv'):
                        final_paths.append(os.path.join(self.parsedir, new_name))
                        #ln trivia: symlink target should be relative to sym link, not where you are.
                        shortcut = os.path.join(self.parsedir, new_name)
                        self.create_symlink(new_name, shortcut)
                
            #delete the output from regular msg parser
            if not self.options.debug:
                shutil.rmtree(expected_msg_parse_dir)
        
        #Clean up
        if not self.options.debug:
            #delete the barcode_splitter output reads files
            subprocess.check_call("rm %s*_read_*" % prefix, shell=True)
            
        log.close()
        
        return final_paths

    def parse(self, use_raw_data_file=None):
        
        #store and return the final relative paths of all the parsed individuals' fastq files
        final_paths = []

        # Convert Illumina reads to format for BWA - Peters program does not create a new file,
        # instead it writes to existing files. So first need to remove the older files.
        # check if output files already exist. If so, erase them.
        if os.path.isfile("./bad_barcodes"): # shorthand-if bad_barcodes file exists, then so do others
            os.remove("./bad_barcodes")
            os.remove("./linker")
            for ind in self.bc:
                fastq_file = 'indiv' + ind[1] + '_' + ind[0]
                os.remove("./" + fastq_file)

        raw_data = use_raw_data_file or self.options.raw_data_file
        stats_file = open(raw_data+'_stats.txt','w')
        
        print "Parsing data (%s) into individual barcode files" % raw_data    
        if self.options.new_parser:
            args = ["python", os.path.join(os.path.dirname(__file__), "grepfqparser.py"),
                raw_data, self.options.barcodes_file, raw_data + '_parsed/']
            print ' '.join(args)
            sys.stdout.flush()
            subprocess.check_call(args)
        else:
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

        total_reads = 0
        for ind in self.bc:
            file_name = 'indiv' + ind[1] + '_' + ind[0]
            fastq_file = raw_data + '_parsed/' + file_name
            number_reads = count_lines(fastq_file)/4
            total_reads += number_reads
            stats_file.write('indiv%s_%s\t%s\n' %(ind[1],ind[0],number_reads))
            #gzip output
            if GZIP_OUTPUT:
                f_in = open(fastq_file, 'rb')
                f_out = gzip.open('%s.gz' % fastq_file, 'wb')
                f_out.writelines(f_in)
                f_out.close()
                f_in.close()
                os.remove(fastq_file)
                fastq_file += '.gz'
            else:
                #TEMP - make sym links to parsed fq files with .fq extensions to they can run in stampy
                #Later versions of stampy will hopefully fix this so we can remove it.
                #ln trivia: symlink target should be relative to sym link, not where you are.
                self.create_symlink(file_name, fastq_file)
            final_paths.append(fastq_file)
        
        stats_file.write('total_reads\t%s' %(total_reads))
        stats_file.close()
    
        self.parsed_time = time.time()
        parsing_time = (self.parsed_time - self.start_time)
        print "Parsing took about %s minutes" %(parsing_time/60)
        return final_paths
    
    def create_symlink(self, file_path, shortcut_path):
        sym_link_args = 'ln -s %s %s.fq' % (file_path, shortcut_path)
        print "create ln with .fq extension:", sym_link_args
        subprocess.check_call([sym_link_args], shell=True)
    
    def _map_w_stampy(self, fastq_file, parent1, parent2, aln_par1_sam, aln_par2_sam, file_par1_log,
        file_par2_log, misc_indiv_log):

        #TEMP: stampy requires .fq extension on our files.  Remove this if stampy fixes it.
        # I'm temporarily creating symlinks in the parse step. (doesn't apply to gz files)
        if not fastq_file.lower().endswith('.gz'):
            fastq_file += '.fq'        
        
        #Align each to sim - output fastq file
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
            #Do the sort (samtools adds .bam suffix to output FYI)
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
            #Do some figuring to get at the fastq file our previous parsing hath made
            fastq_file_name = 'indiv' + ind[1] + '_' + ind[0]
            fastq_file = './' + raw_data + '_parsed/' + fastq_file_name
            #But wait, maybe the parsed fastq file was gzipped so try if it it doesn't exist.
            if not os.path.exists(fastq_file):
                fastq_file += '.gz'
                assert os.path.exists(fastq_file),"File %s could not be found.  Something has gone wrong." % fileq_file

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
            put_back_command = self.options.debug and 'cp' or 'mv' #means 'cp' if DEBUG else 'mv'
            
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
