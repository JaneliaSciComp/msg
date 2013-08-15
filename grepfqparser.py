#!/usr/bin/env python


"""
Parse fastq file containing inline barcodes into separate files

DEPENDENCIES

Python 2.7

USAGE

python grepfqparser.py <input_fastq> <barcode_file> <output_folder>


e.g.

python grepfqparser.py test.fq bc.txt parsed_files


gzipped files will be detected automatically (by .gz suffix) and gunzipped prior to parsing.
This is much faster than parsing on the gzipped files.
The gunzipped file is deleted after parsing.

David L. Stern
Janelia Farm Research Campus
28 May 2013

(There is also a standalone version of this file hosted here: https://github.com/dstern/grepfqparse)

/*
 * Copyright 2013 Howard Hughes Medical Institute.
 * All rights reserved.
 * Use is subject to Janelia Farm Research Campus Software Copyright 1.1
 * license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html ).
 */

"""


import os,sys
import re
import getopt
import subprocess
import uuid

def main():
        """Perform the parsing.  
        Note that this will normally be running in parallel so 
        for any generically named files we write out, we append the first barcode 
        (or a uuid) to them."""
        #parse command line options
        try:
                opts, arg = getopt.getopt(sys.argv[1:],"ht:", ["help"])
        except getopt.error, msg:
                print msg
                print "for help use --help"
                sys.exit(2)
        # process options
        offset = 0
        for o, a in opts:
                if o == '-t':
                    offset = int(a)
                    print "using offset", offset
                if o in ("-h", "--help"):
                        print __doc__
                        sys.exit(0)
        if len(arg) < 3:
                print "\nUsage: python grepfqparser.py [options] <input_fastq> <barcode_file> <output_folder>\m"                
                sys.exit(0)
        #process arguments
        fqFile = arg[0]
        bcFile = arg[1]
        OutFolder = arg[2]
        if len(arg) > 3:
                gzbool = arg[3]
                gzbool = gzbool.upper()
        else:
                gzbool = "YES"


        if os.path.isdir(OutFolder):
                print "Directory %s exists" %(OutFolder)
        else:
                os.mkdir(OutFolder)
        
        print "fastq file = %s" %(fqFile)
        print "barcode file = %s" %(bcFile)
        print "Folder of parsed sequences = %s" %(OutFolder)
        
        #get a unique id to append to file names to avoid overwriting
        #files from other parellel processes.
        uid = str(uuid.uuid1())
        
        y = fqFile.split('.')
        if y[-1] == 'gz':
                gzbool  = "YES"
        else:
                gzbool  = "NO"
                
        print "fq file gzipped? = %s" %(gzbool)
        if gzbool == "YES":
                print "unzipping file"
                tempfq = open("tempfq-"+uid,'w')
                errlog = open("errlog1-"+uid,'w')
                cmd = 'gunzip -c %s' % (fqFile)
                subprocess.check_call(cmd,shell=True,stdout=tempfq,stderr=errlog)
                fqFile = "tempfq-"+uid
                tempfq.close()
                errlog.close()

        bc = open(bcFile,'r')
        bc_for_file_names = ''
        for line in bc:
                lineItems = line.split()
                barcode = lineItems[0]
                barcode_up = barcode.upper()
                #Use first barcode in file for naming files to avoid overwriting files of other processes.
                if not bc_for_file_names:                
                    bc_for_file_names = re.sub(r'[^0-9a-zA-Z]', '', barcode_up)
                name = lineItems[1]
                print barcode
                parsed_file_name = str(OutFolder + "/indiv" + name + "_" + barcode)
                parsed_file = open(parsed_file_name,'w')
                errlog = open("errlog2-" + bc_for_file_names,'w')
                #First grep finds lines starting with the barcode and includes the line above, and two lines below each match
                #Pipe into grep again to filter out -- between matches which some (versions??) of grep insert
                #Pipe into sed to remove barcodes (optional offset) and associated quality scores from each line
                cmd = """grep -B 1 -A 2 ^%s %s | grep -v "^--$" | sed '2~2s/^%s//g'""" % (barcode_up, fqFile, '.'*(len(barcode_up)+offset))
                try:
                    failed = subprocess.call(cmd, shell=True,stdout=parsed_file, stderr=errlog)
                finally:
                    errlog.close()
                    parsed_file.close()
        bc.close()
  
        """Now collect all unparsed reads"""
        bc = open(bcFile,'r')
        allbc = []
        for line in bc:
                lineItems = line.split()
                barcode = lineItems[0]                
                allbc.append(barcode.upper())
        bc.close()
        
        #This code is commented out since we're typically only running on one 
        #barcode at a time and would thus show every other barcode as a non-match
        '''
        """save barcodes in new file -- bconly -- with caret at front, use this as search file"""
        output_bc_file = open("bcOnly-" + bc_for_file_names,'w')
        for line in allbc:
                output_bc_file.write("^" + line + "\n")
        output_bc_file.close()
        bconly  = "bcOnly-" + bc_for_file_names
        nomatch_file = open(OutFolder + "/nomatches-" + bc_for_file_names,'w')
        errlog = open("errlog4-" + bc_for_file_names,'w')
        cmd = "awk 'NR%%4==2' %s | grep -f %s -v" % (fqFile, bconly)
        try:
            nomatch = subprocess.call(cmd, shell=True, stdout=nomatch_file,stderr=errlog)
        finally:
            errlog.close()
            nomatch_file.close()
        '''
        
        """delete tempfq, the gunzipped original file"""
        cmd = 'rm tempfq-%s bcOnly-%s errlog1-%s errlog2-%s' % (
            uid,bc_for_file_names,uid,bc_for_file_names)
        subprocess.call(cmd,shell=True)
        
        
if __name__ == "__main__":
        sys.exit(main())
