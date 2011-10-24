
"""
Make sure each of the important output files contains all of the barcodes.

Expects to be run from directory with barcodes file and output files.
(This sounds bad but all of MSG already expects this.)
"""

import sys
import re

FILES_TO_VALIDATE = ('ancestry-probs-par1.tsv', 'ancestry-probs-par1par2.tsv',
        'ancestry-probs-par2.tsv');

def get_output_file_barcodes(file_path):
    barcodes = set()
    for i,line in enumerate(open(file_path)):
        if i==0:
            continue #skip header line
        #read first column, part after _
        barcodes.add(re.split(r'\s+',line)[0].split('_')[1])
    return barcodes

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

if __name__ == '__main__':
    sys.exit(validate_output_files(sys.argv[1]))

