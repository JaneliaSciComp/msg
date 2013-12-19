
"""
Make sure each of the important output files contains all of the barcodes.

Expects to be run from directory with barcodes file and output files.
(This sounds bad but all of MSG already expects this.)
"""

import sys
import re
from subprocess import Popen, PIPE

FILES_TO_VALIDATE = ('ancestry-probs-par1.tsv', 'ancestry-probs-par1par2.tsv',
        'ancestry-probs-par2.tsv');
#We should see all barcodes when we grep for these patterns in the specified folders
FOLDERS_PATTERNS = {
    'hmm_data':('hmmdata','pileup'),
    'hmm_fit':('hmmprob.pdf','hmmprob.RData'),
}


def get_output_file_barcodes(file_path):
    barcodes = set()
    for i,line in enumerate(open(file_path)):
        if i==0:
            continue #skip header line
        #read first column, part after _
        barcodes.add(re.split(r'\s+',line)[0].split('_')[-1])
    return barcodes

def run_ls_grep(folder,pattern):
    cmd = 'ls -lR %s | grep "%s"' % (folder,pattern)
    #print cmd
    p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = p.communicate()  
    return stdout

def verify_hmm_fit_images():
    stdout = run_ls_grep('hmm_fit_images','pdf')
    sizes = [int(line.split()[4]) for line in stdout.splitlines() if int(line.split()[4])>2000]
    assert len(sizes)>=4,"Missing or too small PDF file in hmm_fit_images"
    stdout = run_ls_grep('hmm_fit_images','bmp')
    sizes = [int(line.split()[4]) for line in stdout.splitlines() if int(line.split()[4])>150000]
    assert len(sizes)>=1,"Missing or too small BMP file in hmm_fit_images"

def verify_all_barcodes_got_important_files(barcodes):
    for folder,patterns in FOLDERS_PATTERNS.items():
        for pattern in patterns:
            #print folder, pattern
            stdout = run_ls_grep(folder, pattern)
            #print stdout
            #print '----'
            for barcode in barcodes:
                if not barcode in stdout:
                    raise ValueError(
                        "A file matching the pattern '%s' for Barcode '%s' not found in any subdirectory of '%s'" % (
                        pattern, barcode, folder))
            for line in stdout.splitlines():
                size = line.split()[4]
                if int(size)<2000:
                    err_msg = "File mentioned in this line is unexpectedly small: " + line
                    sys.stderr.write(err_msg + '\n')

def validate_output_files(barcode_files_path):
    barcodes = set([re.split(r'\s+',line)[0] for line in open(barcode_files_path)])
    #print barcodes
    for output_file_path in FILES_TO_VALIDATE:
        output_barcodes = get_output_file_barcodes(output_file_path)
        #print output_barcodes
        if output_barcodes - barcodes:
            msg = "Validation Error: Output file %s has extra barcodes %s not in barcodes file." % (
                output_file_path, ','.join(output_barcodes - barcodes))
            raise ValueError(msg)
        elif barcodes - output_barcodes:
            msg = "Validation Error: Output file %s is missing barcodes %s." % (
                output_file_path, ','.join(barcodes - output_barcodes))
            raise ValueError(msg)
        else:
            pass
    verify_all_barcodes_got_important_files(barcodes)
    verify_hmm_fit_images()
    
if __name__ == '__main__':
    sys.exit(validate_output_files(sys.argv[1]))

