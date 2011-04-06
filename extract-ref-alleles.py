#!/usr/bin/env python

# [[file:~/Work/simsec/org/simsec.org::*Filter%20reads][block-25]]
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import pysam
from itertools import *
import re, sys, os
from cmdline.cmdline import CommandLineApp
from msglib import *
import copy

class App(CommandLineApp):
    def __init__(self):
        CommandLineApp.__init__(self)
        
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

        op.add_option('--verbosity', action="store_true", dest='verbosity', default=False,
                      help="make lots of noise [default]")

#        op.add_option("-q", "--quiet", action="store_false", dest="verbose",
#                      help="be vewwy quiet (I'm hunting wabbits)")
		  

    def main(self):
        ref_seqs = dict(par1=SeqIO.to_dict(SeqIO.parse(open(self.options.parent1),"fasta")), par2=SeqIO.to_dict(SeqIO.parse(open(self.options.parent2),"fasta")))

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

        print bcolors.OKBLUE + 'Opened SAM files for processing' + bcolors.ENDC;
        sim_refs = dict(zip(sim_samfile.references, [{} for i in range(sim_samfile.nreferences)]))
        sec_refs = dict(zip(sim_samfile.references, [{} for i in range(sim_samfile.nreferences)]))
        refs = dict(par1=sim_refs, par2=sec_refs)
        orths = copy.deepcopy(refs)

        both_samfiles = izip(sim_samfile.fetch(), sec_samfile.fetch())

        species = ['par1','par2']
        refsp = ['par1','par2'] ## species for which reference allele information is being extracted from MD
        self.options.chroms = self.options.chroms.split(',')

        print('Filtering reads and extracting %s reference allele information' % refsp)
        print(sim_samfilename)
        print(sec_samfilename)
        print('Processing contigs: %s' % ' '.join(self.options.chroms))

        i = 0
        read = {}
        for (read['par1'], read['par2']) in both_samfiles:
            if not i % 1e5: print i
            i += 1
  
            ref = sim_samfile.references[read['par1'].rname]
            par2ref = dict(par1=sim_samfile.references[read['par1'].rname],par2=sec_samfile.references[read['par2'].rname])

            if self.options.chroms != ['all'] and ref not in self.options.chroms: continue

            if not (read['par1'].flag in [0,16] and read['par2'].flag in [0,16]):
                if self.options.verbosity: print 'read %d %s flags (parent1: %d, parent2: %d) not in {0,16}: omitting' % (
                    i, read['par1'].qname, read['par1'].flag, read['par2'].flag)
                continue
            if read['par1'].flag == 4 or read['par2'].flag == 4: continue

            seq_forward = dict(zip(species, [Seq(read[sp].seq, IUPAC.ambiguous_dna) for sp in species]))
            if read['par2'].flag != read['par1'].flag:
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

                record_reference_alleles(
							 refs['par1'][ref], refs['par2'][ref], 
							 orths['par1'][ref], orths['par2'][ref], 
							 read['par1'], read['par2'], 
							 seq_forward['par1'], 
							 ref_seqs['par1'][par2ref['par1']].seq, ref_seqs['par2'][par2ref['par2']].seq, 
							 par2ref['par1'], par2ref['par2'], 
							 self.options.verbosity)

#                for sp in refsp:
#                    record_reference_alleles(refs[sp][ref], orths[sp][ref], read['par1'], read[sp], seq_forward[sp], ref_seqs[sp][par2ref[sp]].seq, par2ref[sp], self.options.verbosity)
  
        print '\n'
        
        refsdir = os.path.join(self.options.outdir, 'refs')
        if not os.path.exists(refsdir): os.mkdir(refsdir)
        for sp in refsp:
            outdir = os.path.join(refsdir, sp)
            if not os.path.exists(outdir): os.mkdir(outdir)

            for ref in refs[sp].keys():
                d = refs[sp][ref] 
                if len(d) == 0: continue
                refs[sp][ref] = sorted(zip(d.keys(), d.values()))
                refs[sp][ref] = ["\t".join(map(str, kv)) for kv in refs[sp][ref]]
                alleles_outfile = open(os.path.join(outdir, ref + '-ref.alleles'), 'w')
                alleles_outfile.write('\n'.join(refs[sp][ref]))
                alleles_outfile.write('\n')

                d = orths[sp][ref] 
                orths[sp][ref] = sorted(zip(d.keys(), d.values()))
                orths[sp][ref] = ["\t".join(map(str, kv)) for kv in orths[sp][ref]]
                orths_outfile = open(os.path.join(outdir, ref + '-orths.alleles'), 'w')
                orths_outfile.write('\n'.join(orths[sp][ref]))
                orths_outfile.write('\n')


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

    if verbosity: print "\n\n" + "*" * 100 + "\n" + str(sim_read.qname) + "\n" + "*" * 100 + "\n"
    if verbosity: print "*** SIM_READ ***"
    if verbosity: print sim_read.cigar
    if verbosity: print cigar_ops_sim
    if verbosity: print "*** READ ***"
    if verbosity: print read.cigar
    if verbosity: print cigar_ops
    if verbosity: print "\n*** READ INFO ***\nORIG ref_pos %d" % ref_pos

    if revComp_flag==1: 
        if verbosity: print "reversing order"
        if verbosity: print cigar_ops
        cigar_ops = tuple(reversed(read.cigar))
#        cigar_ops_sim = tuple(reversed(sim_read.cigar))
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
    reconstruct_ref1 = ''
    reconstruct_ref2 = ''

    while index_sim < len(str(sim_ref_seq)):

        while str(sim_ref_seq[index_sim]) == "_" and not str(ref_seq[index_ref]) == "_":
            index_sim += 1
            pos += 1
            reconstruct_ref1 += "."
            reconstruct_ref2 += "."
            if verbosity: print "   \t----- (1 sim %s ref %s) updating index_ref %d and ref_pos %d" % (str(sim_ref_seq[index_sim]),str(ref_seq[index_ref]),index_sim,pos)

        while str(ref_seq[index_ref]) == "_" and not str(sim_ref_seq[index_sim]) == "_":
            index_ref += 1
            if revComp_flag==1: ref_pos -= 1
            else: ref_pos += 1
            reconstruct_ref1 += "."
            reconstruct_ref2 += "."
            if verbosity: print "   \t----- (2 sim %s ref %s) updating index_sim %d and pos %d" % (str(sim_ref_seq[index_sim]),str(ref_seq[index_ref]),index_ref,ref_pos)


        if str(sim_ref_seq[index_sim]) == "_" and str(ref_seq[index_ref]) == "-":
            index_ref += 1
            pos += 1
            if verbosity: print "   \t----- (1*) updating index_ref %d and pos %d" % (index_ref,pos)
        
        elif str(ref_seq[index_ref]) == "_" and str(sim_ref_seq[index_sim]) == "-": 
            index_sim += 1
            if str(ref_seq[index_ref]) != "-":
                if revComp_flag==1: ref_pos -= 1
                else: ref_pos += 1
            if verbosity: print "   \t----- (2*) updating index_sim %d and ref_pos %d" % (index_sim,ref_pos)


        #print "\tsim %d %d: %s\tsec %d %d: %s" % (pos,index_sim,sim_ref_seq[index_sim],ref_pos,index_ref,ref_seq[index_ref])
        #if sim_ref_seq[index_sim] != ref_seq[index_ref]: print "***" 

        if (not str(sim_ref_seq[index_sim]) in ("_","-")) and (not str(ref_seq[index_ref]) in ("_","-")):

            alleles_par1[pos] = par1_ref_seq[pos-1:pos]
            orths_par1[pos] = contig + "\t" + str(pos) + "\t0" + "\t" + str(alleles_par1[pos])

            alleles_par2[pos] = par2_ref_seq[ref_pos-1:ref_pos]
            if revComp_flag: alleles_par2[pos] = par2_ref_seq[ref_pos-1:ref_pos].complement()
            orths_par2[pos] = scaffold + "\t" + str(ref_pos) + "\t" + str(revComp_flag) + "\t" + str(alleles_par2[pos])

            reconstruct_ref1 += alleles_par1[pos]
            reconstruct_ref2 += alleles_par2[pos]

        else: 
            reconstruct_ref1 += "."
            reconstruct_ref2 += "."

        if str(sim_ref_seq[index_sim]) != "-":
            pos += 1

        if str(ref_seq[index_ref]) != "-":
            if revComp_flag==1: ref_pos -= 1
            else: ref_pos += 1

        index_sim += 1
        index_ref += 1

    if verbosity: print "\n*** RECONSTRUCTED ***\n      read  %s\nsim_refseq  %s\n    refseq  %s\n    recons1 %s\n    recons2 %s" % (read_like_sec,sim_ref_seq,ref_seq,reconstruct_ref1,reconstruct_ref2)


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
    try:
        App().run()
    except Exception, e:
        print bcolors.FAIL + 'ERROR in extract-ref-alleles:\n' + bcolors.ENDC
        print '%s' % e
        sys.exit(2)
# block-25 ends here
