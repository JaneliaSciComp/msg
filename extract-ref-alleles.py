#!/usr/bin/env python

# [[file:~/Work/simsec/org/simsec.org::*Filter%20reads][block-25]]
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


class App(CommandLineApp):
    def __init__(self):
        CommandLineApp.__init__(self)
        #self.force_exit = False #enable this when profiling
        
        op = self.option_parser
        op.set_usage('usage: extract_ref_alleles -i indiv -d samdir -o outdir --parent1 parent1_genome --parent2 parent2_genome --chroms csv_list --verbosity')

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

        #I suggest not using this first since it is faster, and using limitmem if
        #you notice you're running out of RAM and hitting swap
        op.add_option('--limitmem', action="store_true", dest='limit_memory', default=False, 
                      help='limit memory usage to size of largest ref/orth file created or 3GB-ish whichever is greater.')

        op.add_option('--verbosity', action="store_true", dest='verbosity', default=False,
                      help="make lots of noise")
	  
    @trace
    def main(self):
        #Clear out refs dir if it exists, create it if it doesn't
        #(It's important to delete any existing files because we open them to write
        #in append mode)
        self.refsdir = os.path.join(self.options.outdir, 'refs')
        if os.path.exists(self.refsdir): 
            shutil.rmtree(self.refsdir)
            assert not os.path.exists(self.refsdir)
        os.mkdir(self.refsdir)

        if not self.options.limit_memory:
            #Disabling garbage collection seemed to take about 7% off of the run time
            #in my testing.
            #Since no large variables are de-referenced or go out of scope anyway
            #there's really no harm to disabling it when not in limit memory mode.
            #(In limit_memory mode we do delete references so we do require 
            #garbage collection.)
            import gc
            gc.disable()

        #Load in species reference genomes (one from each parent.)
        ref_seqs = dict(par1=SeqIO.to_dict(SeqIO.parse(open(self.options.parent1),"fasta")), par2=SeqIO.to_dict(SeqIO.parse(open(self.options.parent2),"fasta")))

        #Manipulates files/paths to get at the barcoded reads for this individual
        #e.g., /home/pinerog/msg_work/MSG_big/SRR071201.fastq_sam_files/aln_indivA1_AAATAG_par1.sam
        sim_samfilename = 'aln_' + self.options.indiv + '_par1.sam'
        sec_samfilename = 'aln_' + self.options.indiv + '_par2.sam'
        #sim_samfilename = 'aln_' + self.options.indiv + '_par1.sam.gz'
        #sec_samfilename = 'aln_' + self.options.indiv + '_par2.sam.gz'
        print(sim_samfilename)
        print(sim_samfilename)

        try:
            sim_samfile = pysam.Samfile(os.path.join(self.options.samdir, sim_samfilename), 'r')
            sec_samfile = pysam.Samfile(os.path.join(self.options.samdir, sec_samfilename), 'r')
        except:
            sim_samfilename += '.gz'
            sec_samfilename += '.gz'
            sim_samfile = pysam.Samfile(os.path.join(self.options.samdir, sim_samfilename), 'r')
            sec_samfile = pysam.Samfile(os.path.join(self.options.samdir, sec_samfilename), 'r')

        sim_outfilename = 'aln_' + self.options.indiv + '_par1-filtered.sam'
        sec_outfilename = 'aln_' + self.options.indiv + '_par2-filtered.sam'


        if not os.path.exists(self.options.outdir): os.mkdir(self.options.outdir)
        if not os.path.exists(self.options.outdir):
            print bcolors.FAIL + 'ERROR in extract-ref-alleles: %s does not exist!' % self.refsdir + bcolors.ENDC
            sys.exit(2)

        sim_outfile = pysam.Samfile(os.path.join(self.options.outdir, sim_outfilename),
                                    'w', template=sim_samfile)
        sec_outfile = pysam.Samfile(os.path.join(self.options.outdir, sec_outfilename),
                                    'w', template=sec_samfile)

        #initialize xxx_refs data structures to an empty dict
        print bcolors.OKBLUE + 'Opened SAM files for processing' + bcolors.ENDC;
        #set up sim ref and sec ref containers - use an empty dict when a key is missing
        #refs, are orths are the only large things kept in memory
        refs = dict(par1=defaultdict(dict), par2=defaultdict(dict))
        orths = copy.deepcopy(refs)

        #Warning: this only iterates to the shortest of the combined lists
        both_samfiles = izip(sim_samfile.fetch(), sec_samfile.fetch())

        species = ['par1','par2']
        refsp = ['par1','par2'] ## species for which reference allele information is being extracted from MD
        self.options.chroms = self.options.chroms.split(',')

        #Optimization: Store these values once instead of re-calc each use
        chroms = self.options.chroms
        verbosity = self.options.verbosity
        omit_flags = set([0,16])
        limit_memory = self.options.limit_memory

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
            par1_rname = read['par1'].rname
            par2_rname = read['par2'].rname
            par1_flag = read['par1'].flag
            par2_flag = read['par2'].flag            

            ref = sim_samfile.references[par1_rname]
            par2ref = dict(par1=sim_samfile.references[par1_rname],par2=sec_samfile.references[par2_rname])

            if chroms != ['all'] and ref not in chroms: continue

            if not (par1_flag in omit_flags and par2_flag in omit_flags):
                if verbosity: print 'read %d %s flags (parent1: %d, parent2: %d) not in {0,16}: omitting' % (
                    i, read['par1'].qname, par1_flag, par2_flag)
                continue
            if par1_flag == 4 or par2_flag == 4: continue

            seq_forward = dict(zip(species, [Seq(read[sp].seq, IUPAC.ambiguous_dna) for sp in species]))
            if par2_flag != par1_flag:
                seq_forward['par2'] = seq_forward['par2'].reverse_complement()
        
            if str(seq_forward[sp]) != str(seq_forward[sp]):
                print read.qname
                print seq_forward['par1']
                print seq_forward['par2']
                raise Exception('Reads on line %d differ' % i)

            try:
                one_best_match = read['par1'].opt('X0') == 1 and read['par2'].opt('X0') == 1
                no_subpoptimal_matches = read['par1'].opt('X1') == 0 and read['par2'].opt('X1') == 0
                mds_same          = read['par1'].opt('MD') == read['par2'].opt('MD')
                cigars_same       = read['par1'].cigar == read['par2'].cigar
                both_species_same = read['par1'].opt('NM') == 0 and read['par2'].opt('NM') == 0
                if both_species_same: assert(mds_same and cigars_same)
                indels            = any([read[sp].opt('XO') > 0 or read[sp].opt('XG') > 0 
                                         for sp in species])
        
                ok = one_best_match and no_subpoptimal_matches
                #ok = one_best_match and no_subpoptimal_matches and not indels
                     ## and not both_species_same

            except: ## KeyError
                ok = False

            if ok:
                sim_outfile.write(read['par1'])
                sec_outfile.write(read['par2'])

                #Fill in the refs and orths dicts
                record_reference_alleles(
							 refs['par1'][ref], refs['par2'][ref], 
							 orths['par1'][ref], orths['par2'][ref], 
							 read['par1'], read['par2'], 
							 seq_forward['par1'], 
							 ref_seqs['par1'][par2ref['par1']].seq, ref_seqs['par2'][par2ref['par2']].seq, 
							 par2ref['par1'], par2ref['par2'], 
							 verbosity)
  
            if (not i % 1e4) and limit_memory:
                #When running in low memory mode, take a break every 1e4 references
                #and offload references calculated so far to the output files
                #so we can clear them from memory.  We'll sort and dedupe them 
                #at the end.
                self.store_alleles_orths(refsp, refs, orths, False)

        print '\n'
        if limit_memory:
            #Write out any remaining refs/orths
            self.store_alleles_orths(refsp, refs, orths, False)
            print "starting sort"
            #load back in the files and finalize them (sort/dedupe)
            self.sort_and_dedupe_ref_orth_files(refsp)
        else:
            #Just store the files once at the end of the run.
            #The will be finalized before being written
            self.store_alleles_orths(refsp, refs, orths, True)

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
        #The only issue is that this program could already
        #be called by qsub 400+ times at once, so if each instance if also calling 
        #100's of qsubs, that might be a problem?
        for sp in refsp:
            target_dir = os.path.join(self.refsdir, sp)
            for fn in os.listdir(target_dir):
                fp = os.path.join(target_dir,fn)
                sort_unique(fp)

        '''
        Misc optimization attempts (included so no-one wastes time repeating these):

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
    def store_alleles_orths(self, refsp, refs, orths, do_finalize_files):
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

def record_reference_alleles(alleles_par1, alleles_par2, orths_par1, orths_par2, sim_read, read, seq_forward, par1_ref_seq, par2_ref_seq, contig, scaffold, verbosity):
    '''Store reference alleles implied by read's MD field. ref_alleles
    is a dict in which the reference alleles are stored. sim_read is
    an AlignedRead object for the same read but aligned against the
    reference genome for which the contig and position information are
    being taken.'''

    pos = sim_read.pos + 1
    ref_pos = read.pos + 1

    revComp_flag = 0
    if sim_read.flag != read.flag and (read.flag == 16 or sim_read.flag == 16): 
        revComp_flag = 1
        ref_pos = read.pos + len(read.seq)

    cigar_ops = read.cigar
    cigar_ops_sim = sim_read.cigar

    if verbosity: 
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
    sim_ref_seq,read_like_sim = updateRefRead(sim_read.qname, cigar_ops_sim, seq_forward, 0, verbosity, "+")

    if verbosity: print "\n*** UPDATING SEC READ ***\n" + "#####"*10
    ref_seq,read_like_sec = updateRefRead(read.qname, cigar_ops, seq_forward, revComp_flag, verbosity, "*")

    if verbosity: print "\n_____FINAL COMPARISON_____\nread     %s\nsim      %s\nlike_sim %s\nsec      %s\nlike_sec %s" % (seq_forward,sim_ref_seq,read_like_sim,ref_seq,read_like_sec)
    if verbosity and len(sim_ref_seq) != len(ref_seq): print "\n_____CHECK ME_____"

    ##################################################
    ### DETERMINE REF ALLELE AND ORTHOLOG
    deletions = re.compile(r"-")
    end_pos = sim_read.pos + len(sim_ref_seq) - len(deletions.findall(str(sim_ref_seq)))
    end_ref_pos = read.pos + len(ref_seq) - len(deletions.findall(str(ref_seq)))

    if verbosity: print("\n*** DETERMINE DIFFS ***\nrevComp_flag %s\npar1_pos %d / read_pos %d - %d\nread   %s\nsimref %s\nrefseq %s\nsimref %s\nsecref %s\nsec_rc %s\n" % 
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
