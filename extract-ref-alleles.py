#!/usr/bin/env python

# [[file:~/Work/simsec/org/simsec.org::*Filter%20reads][block-25]]
import gc
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import pysam
from itertools import *
import re, sys, os
from collections import defaultdict
from cmdline.cmdline import CommandLineApp
from msglib import *
import copy
import shutil
import gzip

##### INTERNAL OPTIONS (for developers) ######

# When true this causes program to periodically offload parts of the data structure
# to disk to keep memory usage lower.  This actually runs slightly faster too.
CONSERVE_MEMORY = True
# How often to offload data structure to disk.  (ignored is CONSERVE_MEMORY is False)
OFFLOAD_MEM_EVERY_N_READS = 5000

##############################################

class App(CommandLineApp):
    def __init__(self):
        CommandLineApp.__init__(self)
        #self.force_exit = False #enable this when profiling
        
        op = self.option_parser
        op.set_usage('usage: extract_ref_alleles -i indiv -d samdir -o outdir --parent1 parent1_genome --parent2 parent2_genome --chroms csv_list --repeat_threshold AS-XS --verbosity')

        op.add_option('-i', '--individual', dest='indiv', type='string', default=None, 
                      help='Name of individual. E.g. for sam file aln_indivA6_AACGAG_sec.sam.gz, indiv would be "indivA6_AACGAG"')

        op.add_option('-d', '--samdir', dest='samdir', type='string', default='', 
                      help='Directory containing SAM files')

        op.add_option('-o', '--outdir', dest='outdir', type='string', default='', 
                      help='Output directory')

        op.add_option('--parent1', dest='parent1', type='string', default=None, 
                      help='Parent 1 genome file to align against')

        op.add_option('--parent2', dest='parent2', type='string', default=None, 
                      help='Parent 2 genome file to align against')

        op.add_option('--chroms', dest='chroms', type='string', default=None, 
                      help='chromosomes desired')

        op.add_option('--verbosity', action="store_true", dest='verbosity', default=False,
                      help="make lots of noise")

        op.add_option('--bwa_alg', dest='bwa_alg', type='string', default='', 
                      help='BWA algorithm used to generate SAM files. It should be "aln" or "bwasw".')

        op.add_option('--use_stampy', dest='use_stampy', type='int', default=0, 
            help='Was STAMPY used to generate SAM files? Set this to 1.')
            
        op.add_option('--repeat_threshold', dest='AS_XS_threshold', type='int', default=6,
            help='Minimum difference between top alignment and suboptimal alignment scores to keep')

    def open_as_pysam(self, filepath):
        try: #Newer pysam deprecates Samfile object in favour of AlignmentFile object
            return pysam.Samfile(filepath, 'r')
        except IOError:
            return pysam.Samfile(filepath + '.gz','r')
        except AttributeError:
            try:
                return pysam.AlignmentFile(filepath, 'r')
            except IOError:
                return pysam.AlignmentFile(filepath + '.gz', 'r')

    @trace
    def remove_non_matching_reads(self, sim_samfilepath, sec_samfilepath, delete_original_files=False):
        """Only keep reads that are in both files (by qname) and don't appear more than 
        once in either file."""   
        sim_removed_count, sec_removed_count = 0,0
        sim_reads, sec_reads = dict(), dict()

        sim_samfile_in = self.open_as_pysam(sim_samfilepath)

        for read in sim_samfile_in.fetch(until_eof=True):  #Newer pysam requires until_eof=True for SAM files, older ignores it
            try: #Newer pysam deprecates .qname in favour of .query_name
                sim_reads[read.query_name] = sim_reads.setdefault(read.query_name,0) + 1
            except:
                sim_reads[read.qname] = sim_reads.setdefault(read.qname,0) + 1
                
        #Close and later re-open input files to save memory and because I don't
        #seem to have access to the seek method in the pysam files to reset the fetch.
        sim_samfile_in.close()
        del sim_samfile_in
        gc.collect()
                                
        sec_samfile_in = self.open_as_pysam(sec_samfilepath)
        for read in sec_samfile_in.fetch(until_eof=True): #Newer pysam requires until_eof=True for SAM files, older ignores it
            try: #Newer pysam deprecates .qname in favour of .query_name
                sec_reads[read.query_name] = sec_reads.setdefault(read.query_name,0) + 1
            except:
                sec_reads[read.qname] = sec_reads.setdefault(read.qname,0) + 1
                
        sec_samfile_in.close()
        del sec_samfile_in
        gc.collect()
            
        both = set(sim_reads.keys()).intersection(set(sec_reads.keys()))
        print "Found %s reads in sim file" % len(sim_reads)
        print "Found %s reads in sec file" % len(sec_reads)
        print "Found %s reads common to both files" % len(both)
        
        #write out only those in both
        
        sim_samfile_in = self.open_as_pysam(sim_samfilepath)
        sim_orig_header = sim_samfile_in.header
        out_sim_samfilepath = sim_samfilepath + '.noncommon.reads.removed.sam'
        try: #Newer pysam deprecates Samfile object in favour of AlignmentFile object
            sim_outfile = pysam.Samfile(out_sim_samfilepath, 'wh', template=sim_samfile_in,
                header=sim_orig_header)
        except AttributeError:
            sim_outfile = pysam.AlignmentFile(out_sim_samfilepath, 'wh', template=sim_samfile_in,
                header=sim_orig_header)
        
        for read in sim_samfile_in.fetch(until_eof=True): #Newer pysam requires until_eof=True for SAM files, older ignores it
            try: #Newer pysam deprecates .qname in favour of .query_name
                if read.query_name in both and (sim_reads[read.query_name] == 1 and sec_reads[read.query_name] == 1):
                    sim_outfile.write(read)
                else:
                    sim_removed_count += 1
            except AttributeError:
                if read.qname in both and (sim_reads[read.qname] == 1 and sec_reads[read.qname] == 1):
                    sim_outfile.write(read)
                else:
                    sim_removed_count +=1
            
        sim_samfile_in.close()
        del sim_samfile_in
        sim_outfile.close()
        del sim_outfile
        gc.collect()
        
        sec_samfile_in = self.open_as_pysam(sec_samfilepath)
        sec_orig_header = sec_samfile_in.header
        out_sec_samfilepath = sec_samfilepath + '.noncommon.reads.removed.sam'
        try: #Newer pysam deprecates Samfile object in favour of AlignmentFile object
            sec_outfile = pysam.AlignmentFile(out_sec_samfilepath, 'wh', template=sec_samfile_in,
                header=sec_orig_header)
        except AttributeError:
            sec_outfile = pysam.Samfile(out_sec_samfilepath, 'wh', template=sec_samfile_in,
                header=sec_orig_header)
        for read in sec_samfile_in.fetch(until_eof=True): #Newer pysam requires until_eof=True for SAM files, older ignores it
            try: #Newer pysam deprecates .qname in favour of .query_name
                if read.query_name in both and (sim_reads[read.query_name] == 1 and sec_reads[read.query_name] == 1):
                    sec_outfile.write(read)
                else:
                    sec_removed_count += 1
            except AttributeError:
                if read.qname in both and (sim_reads[read.qname] == 1 and sec_reads[read.qname] == 1):
                    sec_outfile.write(read)
                else:
                    sec_removed_count +=1 
            
        sec_samfile_in.close()
        del sec_samfile_in
        sec_outfile.close()
        del sec_outfile
        del both
        del sec_reads
        del sim_reads
        gc.collect()
        
        print "Removed %s par1 reads and %s par2 reads" % (sim_removed_count, sec_removed_count)
        
        #Don't delete original since it would prevent us from re-running
        #if delete_original_files:
        #    os.remove(sim_samfilepath)
        #    os.remove(sec_samfilepath)
        
        return out_sim_samfilepath, out_sec_samfilepath
    
    @trace
    def main(self):

        print "\n",get_free_memory(), "mB of free memory at program start"
        print "Conserving Memory mode: %s" % str(CONSERVE_MEMORY)

        #Clear out refs dir if it exists, create it if it doesn't
        #(It's important to delete any existing files because we open them to write
        #in append mode)
        self.refsdir = os.path.join(self.options.outdir, 'refs')
        if os.path.exists(self.refsdir): 
            shutil.rmtree(self.refsdir)
            assert not os.path.exists(self.refsdir)
        os.mkdir(self.refsdir)

        #Disabling garbage collection seemed to take about 7% off of the run time
        #in my testing.
        gc.disable()

        #Load in species reference genomes (one from each parent.)
        def open_file(path):
            if os.path.splitext(path)[1] in ('.gz','.gzip'):
                return gzip.open(path)
            else:
                return open(path)
        ref_seqs = dict(par1=SeqIO.to_dict(SeqIO.parse(open_file(self.options.parent1),"fasta")), 
            par2=SeqIO.to_dict(SeqIO.parse(open_file(self.options.parent2),"fasta")))

        #Manipulates files/paths to get at the barcoded reads for this individual
        #e.g., /home/pinerog/msg_work/MSG_big/SRR071201.fastq_sam_files/aln_indivA1_AAATAG_par1.sam
        sim_samfilename = 'aln_' + self.options.indiv + '_par1.sam'
        sec_samfilename = 'aln_' + self.options.indiv + '_par2.sam'
        sim_samfilepath = os.path.join(self.options.samdir, sim_samfilename)
        sec_samfilepath = os.path.join(self.options.samdir, sec_samfilename)
        
        #sim_samfilename = 'aln_' + self.options.indiv + '_par1.sam.gz'
        #sec_samfilename = 'aln_' + self.options.indiv + '_par2.sam.gz'
        print(sim_samfilename)
        print(sec_samfilename)

        #bwasw throws away bad reads so we can end up with different numbers in each read.
        #we run this to remove any reads that are not in both.
        # (Also use for STAMPY generated SAM files since the reads there also seem to differ.
        # You will also need to run this if you quality trim the alignments.
        # I simply always run it, but if you're running bwa aln, and don't quality trim the alignments
        # you could save some time by not running it.
        sim_samfilepath, sec_samfilepath = self.remove_non_matching_reads(sim_samfilepath, sec_samfilepath,
            True)
        gc.collect()
        
        sim_samfile = self.open_as_pysam(sim_samfilepath)
        sec_samfile = self.open_as_pysam(sec_samfilepath)

        sim_outfilename = 'aln_' + self.options.indiv + '_par1-filtered.sam'
        sec_outfilename = 'aln_' + self.options.indiv + '_par2-filtered.sam'


        if not os.path.exists(self.options.outdir): os.mkdir(self.options.outdir)
        if not os.path.exists(self.options.outdir):
            print bcolors.FAIL + 'ERROR in extract-ref-alleles: %s does not exist!' % self.refsdir + bcolors.ENDC
            sys.exit(2)

        try: #Newer pysam deprecates Samfile object in favour of AlignmentFile object
            sim_outfile = pysam.Samfile(os.path.join(self.options.outdir, sim_outfilename),
                                        'w', template=sim_samfile)
            sec_outfile = pysam.Samfile(os.path.join(self.options.outdir, sec_outfilename),
                                        'w', template=sec_samfile)
        except:
            sim_outfile = pysam.AlignmentFile(os.path.join(self.options.outdir, sim_outfilename),
                                        'w', template=sim_samfile)
            sec_outfile = pysam.AlignmentFile(os.path.join(self.options.outdir, sec_outfilename),
                                        'w', template=sec_samfile)

        #initialize xxx_refs data structures to an empty dict
        print bcolors.OKBLUE + 'Opened SAM files for processing' + bcolors.ENDC;
        #set up sim ref and sec ref containers - use an empty dict when a key is missing
        #refs, are orths are the only large things kept in memory
        refs = dict(par1=defaultdict(dict), par2=defaultdict(dict))
        orths = copy.deepcopy(refs)

        #Warning: this only iterates to the shortest of the combined lists
        both_samfiles = izip(sim_samfile.fetch(until_eof=True), sec_samfile.fetch(until_eof=True)) #Newer pysam requires until_eof=True for SAM files, older ignores it

        species = ['par1','par2']
        refsp = ['par1','par2'] ## species for which reference allele information is being extracted from MD
        self.options.chroms = self.options.chroms.split(',')

        #Optimization: Store these values once instead of re-calc each use
        chroms = self.options.chroms
        verbosity = self.options.verbosity
        omit_flags = set([0,16]) #I believe these are actually "don't omit" flags
        need_sort_and_dedupe = False
        bwa_alg = self.options.bwa_alg.lower()
        use_stampy = self.options.use_stampy

        print('Filtering reads and extracting %s reference allele information' % refsp)
        print(sim_samfilename)
        print(sec_samfilename)
        print('Processing contigs: %s' % ' '.join(self.options.chroms))

        i = 0
        read = {}
        #Go through each reference, keeping only those that match and are valid
        for (read['par1'], read['par2']) in both_samfiles:
            if not i % 1e5: print i
            i += 1
            #Optimization: Store these values once instead of re-calc each use
            try: #Newer pysam deprecates .rname in favour of .reference_id
                par1_rname = read['par1'].reference_id
                par2_rname = read['par2'].reference_id
            except:
                par1_rname = read['par1'].rname
                par2_rname = read['par2'].rname
            par1_flag = read['par1'].flag
            par2_flag = read['par2'].flag            

            ref = sim_samfile.references[par1_rname]
            par2ref = dict(par1=sim_samfile.references[par1_rname],par2=sec_samfile.references[par2_rname])

            if chroms != ['all'] and ref not in chroms: continue

            if not (par1_flag in omit_flags and par2_flag in omit_flags):
                if verbosity:
                    try: #Newer pysam deprecates .qname in favour of .query_name
                        print 'read %d %s flags (parent1: %d, parent2: %d) not in {0,16}: omitting' % (
                        i, read['par1'].query_name, par1_flag, par2_flag)
                    except:
                        print 'read %d %s flags (parent1: %d, parent2: %d) not in {0,16}: omitting' % (
                        i, read['par1'].qname, par1_flag, par2_flag)
                
                continue
            if par1_flag == 4 or par2_flag == 4: continue

            try: #Newer pysam deprecates .seq in favour of .query_sequence
                seq_forward = dict(zip(species, [Seq(read[sp].query_sequence, IUPAC.ambiguous_dna) for sp in species]))
            except:
                seq_forward = dict(zip(species, [Seq(read[sp].seq, IUPAC.ambiguous_dna) for sp in species]))
            if par2_flag != par1_flag:
                seq_forward['par2'] = seq_forward['par2'].reverse_complement()

            if str(seq_forward['par1']) != str(seq_forward['par2']):
                #Output the read that has identical name but different sequence from itself between parental files:
                try:
                    print read['par1'].query_name
                    print read['par2'].query_name
                except:
                    print read['par1'].qname
                    print read['par2'].qname
                print seq_forward['par1']
                print seq_forward['par2']
                raise Exception('Reads on line %d differ' % i)

            #Check for valid reads.
            if bwa_alg == 'bwasw' or use_stampy == 1:
                #For SAM files generated from BWASW or STAMPY algorithms we don't get the same tags.  
                #For now let reads through regardless of one best match or no suboptimal matches
                ok = True
            elif bwa_alg == 'mem':
                #SAM files generated by BWA-MEM do not get the same X* tags as produced by BWA-aln.
                #We can use the AS and XS tags to check for repetitive mapping, and the CIGAR string
                # to check for indels.
                #AS = Alignment Score for reported alignment
                #XS = alignment score for Suboptimal alignments
                #MAPQ=0 also implies multiple mapping for BWA, but only if scores are same
                #i.e. if scores of multiple hits are diff., MAPQ is small, but not 0
                #For divergent parents, we don't bother filtering on identical MDs, CIGARS, or NMs
                #since the read will likely have different numbers of mismatches with each parent,
                #hence different NM, MD, and perhaps CIGAR.
                
                #If CIGAR contains I or D characters, there is an indel
                cigar_indel_chars = ['I', 'D']
                cigar_indel_op_indices = [1, 2]
                try: #Newer pysam provides CIGAR string as .cigarstring
                    indels = cigar_indel_chars in read['par1'].cigarstring or cigar_indel_chars in read['par2'].cigarstring
                except:
                    indels = len([tuple for tuple in read['par1'].cigar if tuple[0] in cigar_indel_op_indices]) > 0 or \
                        len([tuple for tuple in read['par2'].cigar if tuple[0] in cigar_indel_op_indices]) > 0
                
                AS_XS_threshold = self.options.AS_XS_threshold #Default of 6 empirically estimated from Dsim data?
                #AS is Alignment Score of primary alignment
                #XS is alignment Score of suboptimal or unchosen alternate alignment
                #If these two values are too close, it indicates the locus is probably a repeat
                try: #If pysam is recent, opt method is deprecated, use get_tag instead
                    par1_AS_XS_diff = read['par1'].get_tag('AS')-read['par1'].get_tag('XS')
                    par2_AS_XS_diff = read['par2'].get_tag('AS')-read['par2'].get_tag('XS')
                except:
                    par1_AS_XS_diff = read['par1'].opt('AS')-read['par1'].opt('XS')
                    par2_AS_XS_diff = read['par2'].opt('AS')-read['par2'].opt('XS')
                
                non_repetitive = par1_AS_XS_diff > AS_XS_threshold and par2_AS_XS_diff > AS_XS_threshold
                #Only allow reads that lack indels and do not map to repetitive regions:
                ok = not indels and non_repetitive
                
            else:
                try: #Try with newer pysam method .get_tag() and attribute .cigartuples
                    one_best_match = read['par1'].get_tag('X0') == 1 and read['par2'].get_tag('X0') == 1
                    no_suboptimal_matches = False
                    try:
                        no_suboptimal_matches = read['par1'].get_tag('X1') == 0 and read['par2'].get_tag('X1') == 0
                    except:
                        no_suboptimal_matches = read['par1'].get_tag('XT') == 'U' and read['par2'].get_tag('XT') == 'U'
                    mds_same = read['par1'].get_tag('MD') == read['par2'].get_tag('MD')
                    cigars_same = read['par1'].cigartuples == read['par2'].cigartuples #Note: Newer pysam deprecates .cigar in favour of .cigartuples
                    both_species_same = read['par1'].get_tag('NM') == 0 and read['par2'].get_tag('NM') == 0
                    if both_species_same:
                        assert(mds_same and cigars_same), "mds and cigars not the same"
                    indels = any([read[sp].get_tag('XO') > 0 or read[sp].get_tag('XG') > 0 for sp in species])
                    ok = one_best_match and no_suboptimal_matches

                except (KeyError, AssertionError), e:
                    print 'Possible tag missing: %s %s' % (read['par1'].qname, e)
                    ok = False
                except: #Account for newer pysam, where opt is deprecated and replaced by get_tag
                    one_best_match = read['par1'].opt('X0') == 1 and read['par2'].opt('X0') == 1
                    no_suboptimal_matches = False
                    try: no_suboptimal_matches = read['par1'].opt('X1') == 0 and read['par2'].opt('X1') == 0
                    except: no_suboptimal_matches = read['par1'].opt('XT')=='U' and read['par2'].opt('XT')=='U'                    
                    mds_same          = read['par1'].opt('MD') == read['par2'].opt('MD')
                    cigars_same       = read['par1'].cigar == read['par2'].cigar
                    both_species_same = read['par1'].opt('NM') == 0 and read['par2'].opt('NM') == 0
                    if both_species_same: assert(mds_same and cigars_same), "mds and cigars not the same"
                    indels            = any([read[sp].opt('XO') > 0 or read[sp].opt('XG') > 0 for sp in species])
                    ok = one_best_match and no_suboptimal_matches

            if ok:
                sim_outfile.write(read['par1'])
                sec_outfile.write(read['par2'])
                
                #Fill in the refs and orths dicts
                try: #Newer pysam deprecates .seq in favour of .query_sequence
                    record_reference_alleles(
                        refs['par1'][ref], refs['par2'][ref],
                        orths['par1'][ref], orths['par2'][ref],
                        read['par1'], read['par2'],
                        seq_forward['par1'],
                        ref_seqs['par1'][par2ref['par1']].query_sequence, ref_seqs['par2'][par2ref['par2']].query_sequence,
                        par2ref['par1'], par2ref['par2'],
                        verbosity)
                except:
                    record_reference_alleles(
                        refs['par1'][ref], refs['par2'][ref], 
                        orths['par1'][ref], orths['par2'][ref], 
                        read['par1'], read['par2'], 
                        seq_forward['par1'], 
                        ref_seqs['par1'][par2ref['par1']].seq, ref_seqs['par2'][par2ref['par2']].seq, 
                        par2ref['par1'], par2ref['par2'], 
                        verbosity)
                    
  
            if (not i % OFFLOAD_MEM_EVERY_N_READS):
                #print "At ref num %s. free memory is %s mB" % (i, get_free_memory())
                if CONSERVE_MEMORY:
                    #Every 1e4 references check if we are low on memory and if so
                    #offload references calculated so far to the output files
                    #so we can clear them from memory.  We'll sort and dedupe them 
                    #at the end.
                    print "Offloading refs/orths to files to limit further memory usage."
                    print "free memory left: %s mB." % get_free_memory()
                    need_sort_and_dedupe = True
                    self.store_and_remove_alleles_orths(refsp, refs, orths, False)
        print '\n'
        if need_sort_and_dedupe:
            #Write out any remaining refs/orths
            self.store_and_remove_alleles_orths(refsp, refs, orths, False)
            #load back in the files and finalize them (sort/dedupe)
            self.sort_and_dedupe_ref_orth_files(refsp)
        else:
            #Just store the files once at the end of the run.
            #The will be finalized before being written
            self.store_and_remove_alleles_orths(refsp, refs, orths, True)

    @trace
    def sort_and_dedupe_ref_orth_files(self, refsp):
        """
        The files written before calling this function will be out of order, 
        and will have some multiple entries for a given position.  This sorts them and 
        removes any duplicates (keeping the last occurance of a duplicate,
        not the first.  That's how this program has worked
        historically, I'm not sure if it's really neccesary.)
        """
        #If you want a speedup here in the future, each of these sorting jobs
        #could be run by qsub.  They are completely independant.
        #The only issue is that extract-ref-alleles.py could already
        #be called by qsub 400+ times at once, so if each instance if also calling 
        #100's of qsubs, that might be a problem?
        for sp in refsp:
            target_dir = os.path.join(self.refsdir, sp)
            for fn in os.listdir(target_dir):
                fp = os.path.join(target_dir,fn)
                sort_unique(fp)

        '''
        Misc optimization attempts (I'm leaving these in, commented out, so no-one wastes time repeating these):

        #Testing w multiprocessing - used 1GB for each new python process :-( and 
        #seeminly no faster
        from multiprocessing import Pool
        p = Pool(10)
        file_paths = []
        for sp in refsp:
            target_dir = os.path.join(self.refsdir, sp)
            for fn in os.listdir(target_dir):
                fp = os.path.join(target_dir,fn)
                file_paths.append(fp)
        p.map(sort_unique, file_paths)

        #Testing with an in process threadpool http://code.activestate.com/recipes/577187/ (r9)
        #doesn't seem to help.
        pool = ThreadPool(100)
        for sp in refsp:
            target_dir = os.path.join(self.refsdir, sp)
            for fn in os.listdir(target_dir):
                fp = os.path.join(target_dir,fn)
                #pool.add_task(sort_unique, fp)        
        pool.wait_completion()
        '''
    
    @trace
    def store_and_remove_alleles_orths(self, refsp, refs, orths, do_finalize_files):
        """
        Write out the alleles and orths stored in refs and orths.
        We append to the files to handle cases where there is more than one ref/position.
        These files still need to be sorted and deduped.  (see def sort_and_dedupe_ref_orth_files)

        This function can be called before all of the refs and orths have been loaded.
        Which is a handy way to limit memory usage since we delete the refs from the 
        two data strcutures after writing the the files.

        If you set do_finalize_files to True, we assume ALL refs and orths are included
        and will go ahead and sort the data before writing the file.
        """
        for sp in refsp:
            outdir = os.path.join(self.refsdir, sp)
            if not os.path.exists(outdir): os.mkdir(outdir)
            
            for ref in refs[sp].keys():
                d = refs[sp][ref]
                if len(d) == 0: continue
                if do_finalize_files:
                    d_items = sorted(d.items())
                else:
                    d_items = d.items()
                kv_items = ["\t".join(map(str, kv)) for kv in d_items]
                #open files as append so we can keep updating each one as new refs show up 
                #(Doesn't matter when do_finalize_files = True)
                #(just make sure any old files are deleted before we start)
                alleles_outfile = open(os.path.join(outdir, ref + '-ref.alleles'), 'a')
                alleles_outfile.write('\n'.join(kv_items))
                alleles_outfile.write('\n')
                del refs[sp][ref]
                
                d = orths[sp][ref]
                if do_finalize_files:
                    d_items = sorted(d.items())
                else:
                    d_items = d.items()
                kv_items = ["\t".join(map(str, kv)) for kv in d_items]
                orths_outfile = open(os.path.join(outdir, ref + '-orths.alleles'), 'a')
                orths_outfile.write('\n'.join(kv_items))
                orths_outfile.write('\n')
                del orths[sp][ref]
        gc.collect()

def record_reference_alleles(alleles_par1, alleles_par2, orths_par1, orths_par2, sim_read, read, seq_forward, par1_ref_seq, par2_ref_seq, contig, scaffold, verbosity):
    '''Store reference alleles implied by read's MD field. ref_alleles
    is a dict in which the reference alleles are stored. sim_read is
    an AlignedRead object for the same read but aligned against the
    reference genome for which the contig and position information are
    being taken.'''

    try: #Newer pysam deprecates .pos in favour of .reference_start
        pos = sim_read.reference_start + 1
        ref_pos = read.reference_start + 1
    except:
        pos = sim_read.pos + 1
        ref_pos = read.pos + 1

    revComp_flag = 0
    if sim_read.flag != read.flag and (read.flag == 16 or sim_read.flag == 16): 
        revComp_flag = 1
        try: #Newer pysam deprecates .seq in favour of .query_sequence and .pos in favour of .reference_start
            ref_pos = read.reference_start + len(read.query_sequence)
        except:
            ref_pos = read.pos + len(read.seq)

    try: #Newer pysam deprecates .cigar in favour of .cigartuples
        cigar_ops = read.cigartuples
        cigar_ops_sim = sim_read.cigartuples
    except:
        cigar_ops = read.cigar
        cigar_ops_sim = sim_read.cigar

    if verbosity:
        try: #Newer pysam deprecates .qname in favour of .query_name, and .cigar in favour of .cigartuples
            print "\n\n" + "*" * 100 + "\n" + str(sim_read.query_name) + "\n" + "*" * 100 + "\n"
            print "*** SIM_READ ***"
            print sim_read.cigartuples
            print cigar_ops_sim
            print "*** READ ***"
            print read.cigartuples
            print cigar_ops
            print "\n*** READ INFO ***\nORIG ref_pos %d" % ref_pos
        except:
            print "\n\n" + "*" * 100 + "\n" + str(sim_read.qname) + "\n" + "*" * 100 + "\n"
            print "*** SIM_READ ***"
            print sim_read.cigar
            print cigar_ops_sim
            print "*** READ ***"
            print read.cigar
            print cigar_ops
            print "\n*** READ INFO ***\nORIG ref_pos %d" % ref_pos

    if revComp_flag==1: 
        if verbosity: 
            print "reversing order"
            print cigar_ops
        try: #Newer pysam deprecates .cigar in favour of .cigartuples
            cigar_ops = tuple(reversed(read.cigartuples))
            #cigar_ops_sim = tuple(reversed(sim_read.cigartuples))
        except:
            cigar_ops = tuple(reversed(read.cigar))
            #cigar_ops_sim = tuple(reversed(sim_read.cigar))
        if verbosity: print cigar_ops

        for (op, op_len) in cigar_ops:
            if op == 1:		ref_pos -= op_len # insertion to reference
            elif op == 2:	ref_pos += op_len # deletion from reference
        if verbosity: print "NEW  ref_pos %d" % ref_pos

    if verbosity: print "flag " + str(revComp_flag) + "\n" + seq_forward


    ####################################################################################################
    if verbosity: print "\n*** UPDATING SIM READ ***\n" + "#####"*10
    try: #Newer pysam deprecates .qname in favour of .query_name
        sim_ref_seq,read_like_sim = updateRefRead(sim_read.query_name, cigar_ops_sim, seq_forward, 0, verbosity, "+")
    except:
        sim_ref_seq,read_like_sim = updateRefRead(sim_read.qname, cigar_ops_sim, seq_forward, 0, verbosity, "+")

    if verbosity: print "\n*** UPDATING SEC READ ***\n" + "#####"*10
    try: #Newer pysam deprecates .qname in favour of .query_name
        ref_seq,read_like_sec = updateRefRead(read.query_name, cigar_ops, seq_forward, revComp_flag, verbosity, "*")
    except:
        ref_seq,read_like_sec = updateRefRead(read.qname, cigar_ops, seq_forward, revComp_flag, verbosity, "*")

    if verbosity: print "\n_____FINAL COMPARISON_____\nread     %s\nsim      %s\nlike_sim %s\nsec      %s\nlike_sec %s" % (seq_forward,sim_ref_seq,read_like_sim,ref_seq,read_like_sec)
    if verbosity and len(sim_ref_seq) != len(ref_seq): print "\n_____CHECK ME_____"

    ##################################################
    ### DETERMINE REF ALLELE AND ORTHOLOG
    deletions = re.compile(r"-")
    try: #Newer pysam deprecates .pos in favour of .reference_start
        end_pos = sim_read.reference_start + len(sim_ref_seq) - len(deletions.findall(str(sim_ref_seq)))
        end_ref_pos = read.reference_start + len(ref_seq) - len(deletions.findall(str(ref_seq)))
    except:
        end_pos = sim_read.pos + len(sim_ref_seq) - len(deletions.findall(str(sim_ref_seq)))
        end_ref_pos = read.pos + len(ref_seq) - len(deletions.findall(str(ref_seq)))

    if verbosity:
        try: #Newer pysam deprecates .pos in favour of .reference_start
            print("\n*** DETERMINE DIFFS ***\nrevComp_flag %s\npar1_pos %d / read_pos %d - %d\nread   %s\nsimref %s\nrefseq %s\nsimref %s\nsecref %s\nsec_rc %s\n" % 
                (str(revComp_flag),pos,ref_pos,end_ref_pos,seq_forward,sim_ref_seq,ref_seq,
                str(par1_ref_seq[sim_read.reference_start:end_pos]),
                str(par2_ref_seq[read.reference_start:end_ref_pos]),
                str(par2_ref_seq[read.reference_start:end_ref_pos].reverse_complement())))
        except:
            print("\n*** DETERMINE DIFFS ***\nrevComp_flag %s\npar1_pos %d / read_pos %d - %d\nread   %s\nsimref %s\nrefseq %s\nsimref %s\nsecref %s\nsec_rc %s\n" % 
                (str(revComp_flag),pos,ref_pos,end_ref_pos,seq_forward,sim_ref_seq,ref_seq,
                str(par1_ref_seq[sim_read.pos:end_pos]),
                str(par2_ref_seq[read.pos:end_ref_pos]),
                str(par2_ref_seq[read.pos:end_ref_pos].reverse_complement())))

    index_sim = 0
    index_ref = 0
    report_diff = 0
    #Operations to reconstruct_ref1 and reconstruct_ref2 have been commented out for 
    #performance.  If you want to debug you should be able to uncomment them.
    #Some "if verbosity: print..." statements have also been commented out for
    #performance.  If you want to debug you should be able to uncomment them.
    #reconstruct_ref1 = ''
    #reconstruct_ref2 = ''

    #Optimization: Store these values once instead of reconverting on each use
    str_sim_ref_seq = str(sim_ref_seq)
    str_ref_seq = str(ref_seq)
    len_str_sim_ref_seq = len(str_sim_ref_seq)
    while index_sim < len_str_sim_ref_seq:

        while str_sim_ref_seq[index_sim] == "_" and not str_ref_seq[index_ref] == "_":
            index_sim += 1
            pos += 1
            #reconstruct_ref1 += "."
            #reconstruct_ref2 += "."
            #if verbosity: print "   \t----- (1 sim %s ref %s) updating index_ref %d and ref_pos %d" % (str(sim_ref_seq[index_sim]),str(ref_seq[index_ref]),index_sim,pos)

        while str_ref_seq[index_ref] == "_" and not str_sim_ref_seq[index_sim] == "_":
            index_ref += 1
            if revComp_flag==1: ref_pos -= 1
            else: ref_pos += 1
            #reconstruct_ref1 += "."
            #reconstruct_ref2 += "."
            #if verbosity: print "   \t----- (2 sim %s ref %s) updating index_sim %d and pos %d" % (str(sim_ref_seq[index_sim]),str(ref_seq[index_ref]),index_ref,ref_pos)

        if str_sim_ref_seq[index_sim] == "_" and str_ref_seq[index_ref] == "-":
            index_ref += 1
            pos += 1
            #if verbosity: print "   \t----- (1*) updating index_ref %d and pos %d" % (index_ref,pos)
        
        elif str_ref_seq[index_ref] == "_" and str_sim_ref_seq[index_sim] == "-": 
            index_sim += 1
            if str_ref_seq[index_ref] != "-":
                if revComp_flag==1: ref_pos -= 1
                else: ref_pos += 1
            #if verbosity: print "   \t----- (2*) updating index_sim %d and ref_pos %d" % (index_sim,ref_pos)


        #print "\tsim %d %d: %s\tsec %d %d: %s" % (pos,index_sim,sim_ref_seq[index_sim],ref_pos,index_ref,ref_seq[index_ref])
        #if sim_ref_seq[index_sim] != ref_seq[index_ref]: print "***" 

        if (not str_sim_ref_seq[index_sim] in ("_","-")) and (not str_ref_seq[index_ref] in ("_","-")):

            alleles_par1[pos] = par1_ref_seq[pos-1:pos]
            orths_par1[pos] = contig + "\t" + str(pos) + "\t0" + "\t" + str(alleles_par1[pos])

            alleles_par2[pos] = par2_ref_seq[ref_pos-1:ref_pos]
            if revComp_flag: alleles_par2[pos] = par2_ref_seq[ref_pos-1:ref_pos].complement()
            orths_par2[pos] = scaffold + "\t" + str(ref_pos) + "\t" + str(revComp_flag) + "\t" + str(alleles_par2[pos])

            #reconstruct_ref1 += alleles_par1[pos]
            #reconstruct_ref2 += alleles_par2[pos]

        else:
            pass 
            #reconstruct_ref1 += "."
            #reconstruct_ref2 += "."

        if str_sim_ref_seq[index_sim] != "-":
            pos += 1

        if str_ref_seq[index_ref] != "-":
            if revComp_flag==1: ref_pos -= 1
            else: ref_pos += 1

        index_sim += 1
        index_ref += 1

    #if verbosity: print "\n*** RECONSTRUCTED ***\n      read  %s\nsim_refseq  %s\n    refseq  %s\n    recons1 %s\n    recons2 %s" % (read_like_sec,sim_ref_seq,ref_seq,reconstruct_ref1,reconstruct_ref2)


####################################################################################################
### pad reads with insertions and deletions
def updateRefRead(read_id, read_cigar_ops, seq_forward, revComp_flag, verbosity, padding_str):

    updated_seq_forward = ''
    refseq = ''
    i = 0
    j = 0
    for (op, op_len) in read_cigar_ops:
        if verbosity: print "i_%d_op_%d_oplen_%d" % (i,op,op_len)

        if op == 0:	
            refseq += seq_forward[i:i+op_len]
            updated_seq_forward += seq_forward[j:j+op_len]

        elif op == 1:	
            refseq += '-' * op_len
            updated_seq_forward += seq_forward[j:j+op_len]

        elif op == 2:
            refseq += '_' * op_len
            updated_seq_forward += padding_str * op_len

        elif op in (4,5):	
            refseq += 'X' * op_len # ignored but present in reference
            updated_seq_forward += seq_forward[j:j+op_len]

        elif op != 6:
            if verbosity: print "%s unknown op %d" % (read_id,op)
            raise Exception('--> %s <-- Non-standard operation in CIGAR (%d,%d)' % (str(read_id),op,op_len))

        if op != 2: i += op_len
        if op != 2: j += op_len
        if verbosity: print "ref  " + refseq + "\nread " + updated_seq_forward

    return (refseq,updated_seq_forward)


if __name__ == '__main__':
    '''
    #Python Profiling code: (set force_exit to false in __init__ above)
    #Uncomment the code below and comment out the normal call to App.run()
    import hotshot, hotshot.stats
    prof = hotshot.Profile("index.prof")
    prof.runcall(App().run)
    prof.close()
    stats = hotshot.stats.load("index.prof")
    #stats.strip_dirs()
    stats.sort_stats('time', 'calls')
    stats.print_stats()
    print '-------Callers------'
    stats.print_callers()
    print '------Callees-------'
    stats.print_callees() #redundant?
    '''
    try:
        App().run()
    except Exception, e:
        print bcolors.FAIL + 'ERROR in extract-ref-alleles:\n' + bcolors.ENDC
        print '%s' % e
        sys.exit(2)
    #'''
# block-25 ends here
