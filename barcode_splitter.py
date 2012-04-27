#!/usr/bin/env python
# from https://bitbucket.org/lance_parsons/paired_sequence_utils/src/355e838b92d0/barcode_splitter.py

""" Split one or more fastq files based on barcode sequence.
"""
from __future__ import division
import gzip
import optparse
import os
import re
import sys
import mapping_functions
#from fastq_utils import fastq_utils

__version__ = "0.1"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2011, Lance Parsons"
__license__ = "BSD 2-Clause License http://www.opensource.org/licenses/BSD-2-Clause"

UNMATCHED = 'unmatched'

def main (argv=None):
    print "barcode_splitter starting run"
    if argv is None:
        argv = sys.argv
    
    usage = "Usage: %prog [options] fastq_read1 [fastq_read2] [fastq_read3]"
    parser = optparse.OptionParser(usage=usage, version='%prog version ' + globals()['__version__'], description=globals()['__doc__'])
    parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
    
    parser.add_option ('--bcfile', help='Tab delimited file with barcodes and sample ids (REQUIRED)')
    
    
    #When used in MSG (Multiplexed shotgun sequencing) the illumina indexes will be parsed out first
    #and those results will be passed through the MSG parser to find the MSG barcodes.  Since each index
    #will have the same barcodes, the MSG parser renames the results with the index.  
    #This option group creates a new barcodes file for MSG to use downstream from parsing with the updated
    #names.
    msg_group = optparse.OptionGroup(parser, "MSG Options")
    msg_group.add_option ('--make_indexed_msg_barcodes_file',help='If you have a MSG barcodes file, this creates a new version of it with indexes',
        action='store_true', default=False)
    msg_group.add_option ('--msg_barcodes', help="Tab delimited seq/label list of MSG barcodes", default='')
    parser.add_option_group(msg_group)
    
    output_group = optparse.OptionGroup(parser, "Output Options")
    output_group.add_option ('--prefix', default='', help='Prefix for output files')
    output_group.add_option ('--suffix', default='', help='Suffix for output files')
    parser.add_option_group(output_group)
    
    barcode_location_group = optparse.OptionGroup(parser, "Barcode Location")
    barcode_location_group.add_option ('--idxread', default=1, type='int', help='Indicate in which read to search for the barcode sequence (default: 1)')
    barcode_location_group.add_option ('--barcodes_at_end', action='store_true', default=False, help='Barcodes are at the end of the index read (default is at the beginning)')
    parser.add_option_group(barcode_location_group)
    
    matching_group = optparse.OptionGroup(parser, "Matching")
    matching_group.add_option ('--mismatches', default=0, type='int',  help='Number of mismatches allowed in barcode matching')
    parser.add_option_group(matching_group)
    
    input_group = optparse.OptionGroup(parser, "Input format")
    input_group.add_option ('--format', default='fastq', help='Specify format for sequence files (fasta or fastq)')
    input_group.add_option('--gzip', action='store_true', default=False, help='Force gzip format for input and output (default is auto based on input file extension)')
    parser.add_option_group(input_group)
    
    
    try:
        (options, args) = parser.parse_args(argv[1:])
        if (not options.make_indexed_msg_barcodes_file) and len(args) < 1:
            parser.error('Must specify at least one sequence file')
        if not options.bcfile:
            parser.error('Must specify a barcodes file with "--bcfile" option')
    except SystemExit: # Prevent exit when calling as function
        return 2
    
    if options.make_indexed_msg_barcodes_file and options.msg_barcodes:
        return make_msg_barcodes_file(options.msg_barcodes, options.bcfile)
    
    # Read barcodes files into dict
    barcode_dict = read_barcodes(options.bcfile)
    total_read_count = 0
    counts = {UNMATCHED: 0}
    for barcode in barcode_dict:
        counts[barcode] = 0
    # TODO Verbose: print barcode_dict
    
    # Determine if we should use gzip for input/output
    basename, extension = os.path.splitext(args[0])
    if extension == '.gz':
        options.gzip = True

    # Open filehandles for each read
    inputs = {}
    outputs = {}
    for i in xrange(0,len(args)):
        if options.gzip:
            inputs[i] = gzip.open(args[i], 'rb')
        else:
            inputs[i] = open(args[i], 'rb')
        first_line = inputs[i].readline().strip()
        id_format = determine_id_format(first_line[1:])
        # TODO print "\nFile %s: %s\nId format: %s" % (i, args[i], id_format)
        inputs[i].seek(0)
    
    # Open filehandles for each barcode (and unmatched), for each read
    for barcode in barcode_dict:
        outputs[barcode] = {}
        for i in inputs:
            if options.gzip:
                outputs[barcode][i] = gzip.open('%s%s_read_%s%s' % (options.prefix, barcode_dict[barcode], i+1, options.suffix), 'wb')
            else:
                outputs[barcode][i] = open('%s%s_read_%s%s' % (options.prefix, barcode_dict[barcode], i+1, options.suffix), 'wb')
    outputs[UNMATCHED] = {}
    for i in inputs:
        if options.gzip:
            outputs[UNMATCHED][i] = gzip.open('%s%s_read_%s%s' % (options.prefix, UNMATCHED, i+1, options.suffix), 'wb')  
        else:
            outputs[UNMATCHED][i] = open('%s%s_read_%s%s' % (options.prefix, UNMATCHED, i+1, options.suffix), 'wb')    
            
    # For each input line in index read, get index sequence
    for index_read in read_fastq(inputs[options.idxread-1]):
        total_read_count += 1
        barcode_length = len(barcode_dict.keys()[0])
        if options.barcodes_at_end:
            index_seq = index_read['seq'][-barcode_length:]
        else:
            index_seq = index_read['seq'][0:barcode_length]
            
        # Get matching barcode(s), if more than one, warn and set to unmatched
        best_match_barcodes = match_barcodes(index_seq, barcode_dict.keys(), options.mismatches)
        if (len(best_match_barcodes) == 1):
            barcode_match = best_match_barcodes[0]
        else:
            barcode_match = UNMATCHED
            if (len(best_match_barcodes) > 1) :
                sys.stderr.write('More than one barcode matches for %s, moving to %s category\n' % (index_read['seq_id'], UNMATCHED))
        counts[barcode_match] += 1
        
        # Get sequence record from each other read, assert id matches
        for readnum in xrange(0,len(inputs)):
            if readnum != options.idxread-1:
                read = read_fastq(inputs[readnum]).next()
                try:
                    assert(match_id(index_read['seq_id'], read['seq_id'], id_format))
                except AssertionError:
                    sys.stderr.write("Id mismatch: %s does not match %s" %(index_read['seq_id'], read['seq_id']))
            else:
                read = index_read
            # Output sequences into barcode/read file    
            outputs[barcode_match][readnum].write(fastq_string(read))
    print "Sample\tBarcode\tCount\tPercent"
    for barcode in sorted(barcode_dict, key=barcode_dict.get):
        print "%s\t%s\t%s\t%.2f%%" % (barcode_dict[barcode], barcode, counts[barcode], (counts[barcode]/total_read_count)*100 )
    print "%s\t%s\t%s\t%.2f%%" % (UNMATCHED, None, counts[UNMATCHED], (counts[UNMATCHED]/total_read_count)*100 )
    return 0

def read_barcodes(filename):
    '''Read barcodes file into dictionary'''
    barcode_dict = {}
    linenum = 0
    filehandle = open(filename, 'rb')
    for line in filehandle:
        linenum += 1
        line = line.strip()
        if line[0] != '#':
            (sample_id, barcode_sequence) = line.split('\t')
            if (sample_id is not None) and (barcode_sequence is not None):
                barcode_dict[barcode_sequence] = sample_id
            else:
                raise Exception("Unable to read barcode from line %s: '%s'" % (linenum, line))
    return barcode_dict

def make_msg_barcodes_file(msg_barcodes, bcfile):
    """
    When used in MSG (Multiplexed shotgun sequencing) the illumina indexes will be parsed out first
    and those results will be passed through the MSG parser to find the MSG barcodes.  Since each index
    will have the same barcodes, the MSG parser renames the results with the index.  
    This function creates a new barcodes file for MSG to use downstream from parsing with the updated
    names.
    """
    new_barcodes = open(msg_barcodes + '.after.index.parsing','w')
    for index_id in read_barcodes(bcfile).values():
        for ind in mapping_functions.read_barcodes(msg_barcodes):
            #For reference, MSG parsing output file_name is composed of: 'indiv' + ind[1] + '_' + ind[0]
            new_ind = ind[:]
            new_ind[1] = new_ind[1] + index_id
            new_barcodes.write('\t'.join(new_ind) + '\n')
    new_barcodes.close()  
    return 0

def match_barcodes(sequence, barcodes, mismatches):
    '''Find closest match(es) in barcodes to specified sequence with max number of mismatches'''
    best_distance = mismatches
    results = []
    for barcode in barcodes:
        if mismatches == 0:
            if (sequence == barcode):
                results.append(barcode)
        else:
            distance = hamming_distance(sequence, barcode)
            if (distance <= best_distance):
                best_distance = distance
                results.append(barcode)
    return results
    
def hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))    

''' Supported types of Fastq IDs '''
ILLUMINA = 'illumina' # CASVA 1.8+, match up to space
STRIPONE = 'stripone' # Illumina CASAVA 1.7 and lower (/1 /2) and NCBI SRA (/f /r), match all but last character
OTHER = 'other'       # Other, match exactly

def determine_id_format(seq_id):
    '''Determine if the id is new illumina, old illumina (/1 /2 ...), sanger (/f /r), or other'''
    
    id_format = None
    # Illumina CASAVA 1.8+ fastq headers use new format
    read_id_regex = re.compile(r'(?P<instrument>[a-zA-Z0-9_-]+):(?P<run_number>[0-9]+):(?P<flowcell_id>[a-zA-Z0-9]+):(?P<lane>[0-9]+):(?P<tile>[0-9]+):(?P<x_pos>[0-9]+):(?P<y_pos>[0-9]+) (?P<read>[0-9]+):(?P<is_filtered>[YN]):(?P<control_number>[0-9]+):(?P<index_sequence>[ACGT]+){0,1}')
    # Old illumina and sanger reads use /1 /2 or /f /r
    strip_one_endings = ['/1', '/2', '/3', '/f', '/r']
    
    if read_id_regex.match(seq_id):
        id_format = ILLUMINA
    elif (seq_id[-2:] in strip_one_endings):
        id_format = STRIPONE
    else:
        id_format = OTHER
    return id_format

def strip_read_from_id(seq_id, id_format=None):
    new_id = seq_id
    if not id_format:
        id_format = determine_id_format(seq_id)
    elif id_format == STRIPONE:
        new_id = seq_id[0:-1]
    elif id_format == ILLUMINA:
        new_id = seq_id.split(' ')[0]
    return new_id

def strip_read_from_id_stripone(seq_id):
    return seq_id[0:-1]

def strip_read_from_id_illumina(seq_id):
    return seq_id.split(' ')[0]

def match_id(id1, id2, id_format=OTHER):
    ''' Return true if id's match using rules for specified format '''
    if id_format == STRIPONE:
        if id1[0:-1] == id2[0:-1]:
            return True
        else:
            return False
    elif id_format == ILLUMINA:
        if (id1.split(' ')[0] == id2.split(' ')[0]):
            return True
        else:
            return False
    elif id1 == id2:
        return True
    else:
        return False

def read_fastq(filehandle):
    ''' Return dictionary with 'seq_id', 'seq', 'qual_id', and 'qual' '''
    record_line = 0
    read_number = 0
    fastq_record = dict()
    for line in filehandle:
        record_line += 1
        if record_line == 1:
            fastq_record['seq_id'] = line.strip()
        elif record_line == 2:
            fastq_record['seq'] = line.strip()
        elif record_line == 3:
            fastq_record['qual_id'] = line.strip()
        elif record_line == 4:
            record_line = 0
            fastq_record['qual'] = line.strip()
            read_number += 1
            yield fastq_record
            

def fastq_string(record):
    return "%s\n%s\n%s\n%s\n" % (record['seq_id'], record['seq'], record['qual_id'], record['qual'])

if __name__ == '__main__':
    sys.exit(main())
