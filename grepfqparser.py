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
import getopt
import subprocess

def main():
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
                print "\nUsage: python grepfqparser.py <input_fastq> <barcode_file> <output_folder>\m"                
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
        
        y = fqFile.split('.')
        if y[-1] == 'gz':
                gzbool  = "YES"
        else:
                gzbool  = "NO"
                
        print "fq file gzipped? = %s" %(gzbool)        
        if gzbool == "YES":
                print "unzipping file"
                tempfq = open("tempfq",'w')
                errlog = open("errlog1",'w')
                cmd = 'gunzip -c %s' % (fqFile)
                subprocess.check_call(cmd,shell=True,stdout=tempfq,stderr=errlog)
                fqFile = "tempfq"
                tempfq.close()
                errlog.close()

        bc = open(bcFile,'r')
        for line in bc:
                lineItems = line.split()
                barcode = lineItems[0]
                barcode_up = barcode.upper()
                name = lineItems[1]
                print barcode
                
                parsed_file_step1_name = str(OutFolder + "/indiv" + name + "_" + barcode + "firstgrep")
                parsed_file_step1 = open(parsed_file_step1_name,'w') 
                errlog = open("errlog2",'w')
                #(note: pipe into sed to remove barcodes and associated quality scores from each line)
                cmd = "grep -B 1 -A 2 ^%s %s | sed '2~2s/^%s//g'" % (barcode_up, fqFile, '.'*(len(barcode_up)+offset))
                failed = subprocess.call(cmd, shell=True,stdout=parsed_file_step1, stderr=errlog)
                errlog.close()
                parsed_file_step1.close()
                if failed:
                    continue
                """grep with -B and -A produces spacer marks '--' in file. Cannot figure out how to suppress these, 
                so remove and paste into new file, then delete original"""
                
                parsed_file_name = str(OutFolder + "/indiv" + name + "_" + barcode)
                parsed_file = open(parsed_file_name,'w')
                errlog = open("errlog3",'w')
                cmd = 'awk "!/^--$/" %s' % (parsed_file_step1_name)
                subprocess.check_call(cmd,shell=True,stdout=parsed_file,stderr=errlog)
                errlog.close()
                parsed_file.close()
                cmd = 'rm %s' % (parsed_file_step1_name)
                subprocess.check_call(cmd,shell=True)
        
        bc.close()
  
        """Now collect all unparsed reads"""
        bc = open(bcFile,'r')
        allbc = []
        for line in bc:
                lineItems = line.split()
                barcode = lineItems[0]                
                allbc.append(barcode.upper())
        bc.close()
        
        """save barcodes in new file -- bconly -- with caret at front, use this as search file"""
        output_bc_file = open("bcOnly",'w')
        for line in allbc:
                output_bc_file.write("^" + line + "\n")
        output_bc_file.close()
        bconly  = "bcOnly"
        nomatch_file = open(OutFolder + "/nomatches",'w')
        errlog = open("errlog4",'w')
        cmd = "awk 'NR%%4==2' %s | grep -f %s -v" % (fqFile, bconly)
        try:
            nomatch = subprocess.call(cmd, shell=True, stdout=nomatch_file,stderr=errlog)
        finally:
	    errlog.close()
            nomatch_file.close()
        
        """delete tempfq, the gunzipped original file"""
        cmd = 'rm tempfq bcOnly errlog1 errlog2 errlog3 errlog4'
        subprocess.call(cmd,shell=True)
        
        
if __name__ == "__main__":
        sys.exit(main())
