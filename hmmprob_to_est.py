
"""
Find all [name]_[barcode]-hmmprob.RData.chrom.[contig].csv files in hmm_fit directory
and create a new CSV with all individuals and contigs.

example file name: indivH12_TTGACG-hmmprob.RData.chrom.2R.csv

File Format Details:
    Creates a new CSV file at out_path with one row for each individual 
    and a column for every single chromosome, and position on the chromosome.
    Columns contain the est value for that indivdual at that position.
    
Usage:
python msg/hmmprob_to_est.py -d hmm_fit -o hmmprob_w_est.csv
"""

import os
import sys
import csv
import glob
import optparse

# -------------- SETTINGS ----------------
# Assumes the files matches this pattern relative to hmm_fit (or other specificied directory)
GLOB_PATTERN = '/*/*-hmmprob.RData.chrom.*.csv'

# ----------------------------------------

def grab_files(dir):
    """Example:
    glob.glob('hmm_fit/*/*-hmmprob.RData.chrom.*.csv')
    ['hmm_fit/indivF11_GTTACG/indivF11_GTTACG-hmmprob.RData.chrom.2R.csv', 'hmm_fit/indivE2_CAGCCG/indivE2_CAGCCG-hmmprob.RData.chrom.2R.csv', 'hmm_fit/indivG7_CTTGCG/indivG7_CTTGCG-hmmprob.RData.chrom.2R.csv', 'hmm_fit/indivH6_GTCCGG/indivH6_GTCCGG-hmmprob.RData.chrom.2R.csv', 'hmm_fit/indivH9_GTCGCG/indivH9_GTCGCG-hmmprob.RData.chrom.2R.csv', 'hmm_fit/indivE12_GTATCG/indivE12_GTATCG-hmmprob.RData.chrom.2R.csv', 'hmm_fit/indivH12_TTGACG/indivH12_TTGACG-hmmprob.RData.chrom.2R.csv', 'hmm_fit/indivA12_AATAAG/indivA12_AATAAG-hmmprob.RData.chrom.2R.csv', 'hmm_fit/indivE8_CTAACG/indivE8_CTAACG-hmmprob.RData.chrom.2R.csv', 'hmm_fit/indivE4_TAGGAG/indivE4_TAGGAG-hmmprob.RData.chrom.2R.csv']
    """
    glob_pattern = dir.strip('/') + GLOB_PATTERN
    return glob.glob(glob_pattern)

def parse_path(path):
    """Get ind name, and chrom from file path"""
    dir, filename = os.path.split(path)
    name_parts = filename.split('.')
    ind_name = name_parts[0].strip('-hmmprob')
    chrom = name_parts[-2]
    return ind_name, chrom

def transform(file_list, out_path):
    """
    Creates a new CSV file at out_path with one row for each individual 
    and a column for every single chromosome, and position on the chromosome.
    Columns contain the est value for that indivdual at that position.
    """
    #First store all the ests by position for each individual in ests_by_ind.
    # { 'ind_name': {('2R',1):est, ('2R',3):est}, 'ind_name2': {...}, }
    ests_by_ind = {}
    
    for path in file_list:
        ind_name, chrom = parse_path(path)
        if not ind_name in ests_by_ind:
            ests_by_ind[ind_name] = {}
        csv_reader = csv.reader(open(path, 'rb'))
        csv_reader.next() #skip header row
        for row in csv_reader:
            pos, count, est = row[1], row[15], row[20]
            ests_by_ind[ind_name][(chrom,pos)] = est

    #Make a sorted list of all chroms/positions
    all_positions = set()
    for dpositions in ests_by_ind.values():
        all_positions |= set(dpositions.keys())
    all_positions = sorted(list(all_positions))
    
    #Write to CSV
    outfile = open(out_path, 'wb')
    outcsv = csv.writer(outfile)
    #Write header row
    outcsv.writerow(['individual'] + ['-'.join(v) for v in all_positions])
    #Write data
    for ind_name, ests_by_chrom_pos in ests_by_ind.items():
        outrow = [ind_name]
        for (chrom, pos) in all_positions:
            outrow.append(ests_by_chrom_pos.get((chrom, pos),'N/A'))
        outcsv.writerow(outrow)
    outfile.close()
    
def main():
    """Parse command line args, and call appropriate functions."""
    usage="""\nusage: %prog [options]\n"""
    parser = optparse.OptionParser(usage=usage)
    #Other option types are int and float, string is default.
    #Note there is also a default parameter.
    parser.add_option('-d','--dir',dest="hmm_fit_dir",type="string")
    parser.add_option('-o','--out',dest="out_path",type="string")
    parser.add_option('-t','--thresh',dest="pna_thresh",type="float",default=.03)
    opts,args=parser.parse_args() #Args taken from sys.argv[1:] by default, parsed using GNU/POSIX syntax.
    if not opts.hmm_fit_dir and opts.out_path:
        parser.error("A directory for locating hmm_fit data and output file path is required.")

    all_files = grab_files(opts.hmm_fit_dir)
    
    transform(all_files, opts.out_path)
    
if __name__=='__main__':
    main()