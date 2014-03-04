
"""
Find all [name]_[barcode]-hmmprob.RData.chrom.[contig].csv files in hmm_fit directory
and create a new CSV with all individuals and contigs.

example file name: indivH12_TTGACG-hmmprob.RData.chrom.2R.csv

(If file name format changes, update parse_path function below)

File Format Details:
    Creates a new CSV file at out_path with one individual per row.
    Columns contain the est value for that indivdual at that position.
    (Columns with less than pnathresh % coverage are omitted.)
    
Usage:
python msg/hmmprob_to_est.py -d hmm_fit -t .03 -o hmmprob_w_est.csv
"""

import os
import sys
import csv
import glob
import optparse

from msglib import trace, get_free_memory

# -------------- SETTINGS ----------------

# Assumes the files matches this pattern relative to hmm_fit (or other specificied directory)
GLOB_PATTERN = '/*/*-hmmprob.RData.chrom.*.csv'

# 0 based column indexes for columns we care about in input CSV files
COL_POS = 1
COL_COUNT = 15 #not used
COL_EST = 20

# ----------------------------------------

def grab_files(dir):
    """Example from Toy data:
    glob.glob('hmm_fit/*/*-hmmprob.RData.chrom.*.csv')
    ['hmm_fit/indivF11_GTTACG/indivF11_GTTACG-hmmprob.RData.chrom.2R.csv', 'hmm_fit/indivE2_CAGCCG/indivE2_CAGCCG-hmmprob.RData.chrom.2R.csv', 
    'hmm_fit/indivG7_CTTGCG/indivG7_CTTGCG-hmmprob.RData.chrom.2R.csv', ...]
    """
    glob_pattern = dir.strip('/') + GLOB_PATTERN
    return glob.glob(glob_pattern)

@trace
def write_csv(d_ests, out_path):
    """ Takes a (filtered) dict (See transform function for explanation of
    what d_ests is and an example of what it could contain.)
    
    Writes it out into a CSV file, putting everything together in one matrix.    
    """
    #Write to CSV
    outfile = open(out_path, 'wb')
    outcsv = csv.writer(outfile)

    #Set up data to fill in, and be written to file
    #Header rows:
    header_row, chrom_row, gen_map_pos_row = ['individual'],[''],['']
    #seed with ind names (Make first column of each data row hold the ind name.)
    csv_data = []
    for d_inds in d_ests.values():
        for ind_name in d_inds:
            csv_data.append((ind_name,))
    #Sort and make sure there are no duplicates
    csv_data = sorted(list(set(csv_data)))
    #change rows from tuples to lists
    csv_data = [list(row) for row in csv_data]
    #print csv_data
    #build an index of our csv_data so we can quickly put an indivudals data in the right place
    r = row_by_ind_name = {}
    for i,row in enumerate(csv_data):
        r[row[0]] = i

    #fill in data
    for chrom, d_inds in d_ests.items():
        #Make a sorted list of all positions in this chromosome
        all_positions = set()
        for d_ests in d_inds.values():
            all_positions |= set(d_ests.keys())
        all_positions = sorted(list(all_positions))
        
        if all_positions: #only include chroms with data
            #Update header rows with these positions / chomosomes
            header_row += ['%s-%s' % (chrom, v) for v in all_positions]
            chrom_row += ([chrom] * len(all_positions))
            gen_map_pos_row += [i+1 for i in range(len(all_positions))]

            #Store actual data to be written
            for ind_name, ests_by_pos in d_inds.items():
                outrow = csv_data[r[ind_name]]
                for pos in all_positions:
                    outrow.append(ests_by_pos.get(pos,'-'))
    
    outcsv.writerow(header_row)
    outcsv.writerow(chrom_row)
    outcsv.writerow(gen_map_pos_row)
    outcsv.writerows(csv_data)
    print "Free memory after writing CSV is %s MB" % get_free_memory()
    outfile.close()

def parse_path(path):
    """Get ind name, and chrom from file path"""
    dir, filename = os.path.split(path)
    name_parts = filename.split('.')
    ind_name = name_parts[0].strip('-hmmprob')
    chrom = name_parts[-2]
    return ind_name, chrom

@trace
def transform(file_list, pnathresh):
    """
    Groups position ests by individual and by chromosome and filters 
    out positions with less than pnathresh % coverage.
    """
    
    #d_ests stores estimates by position by individual by chromosome.
    #example:
    #    {'2R': {'indivA12_AATAAG': {'1000992': '1',
    #                                '10065531': '3',
    #                                '9987712': '1'},
    #            'indivE12_GTATCG': {'10002269': '3',
    #                                '10022498': '3',
    #                                '10079005': '3'},
    #                                },
    #     '3R': ...
    #        }
    d_ests = {}
    chrom_pos_count = {} #count of individuals with a given (chrom,position)
    
    #Fill up data structure from all files
    for path in file_list:
        ind_name, chrom = parse_path(path)
        if not chrom in d_ests:
            d_ests[chrom] = {}
        if not ind_name in d_ests[chrom]:
            d_ests[chrom][ind_name] = {}
        csv_reader = csv.reader(open(path, 'rb'))
        csv_reader.next() #skip header row
        for row in csv_reader:
            pos, count, est = row[COL_POS], row[COL_COUNT], row[COL_EST]
            d_ests[chrom][ind_name][pos] = est
            chrom_pos_count[(chrom,pos)] = chrom_pos_count.get((chrom,pos),0) + 1

    print "(mid transform function) Free memory now is %s MB" % get_free_memory()

    #Remove positions with less individuals than pna thresh %
    #(example: If pna thresh is .1, that means for a given chromosome location
    #we'd throw out the whole position if it exists for less than 10% of individuals)
    num_inds = max([len(d_inds) for d_inds in d_ests.values()])
    print "There are %s individuals" % num_inds
    count_thresh = int(round(pnathresh * num_inds))
    print "Will throw out chrom/positions with less than %s individuals." % count_thresh
    print "(that's int(round(pna_thresh %s * %s individuals)) = %s )" % (pnathresh, num_inds, count_thresh)

    for chrom, d_inds in d_ests.items():
        for ind_name, ests_by_pos in d_inds.items():
            for pos in ests_by_pos.keys():
                if chrom_pos_count[(chrom,pos)] < count_thresh:
                    del d_ests[chrom][ind_name][pos]                
    return d_ests

@trace
def main():
    """Parse command line args, and call appropriate functions."""
    usage="""\nusage: %prog [options]\n"""
    parser = optparse.OptionParser(usage=usage)
    #Other option types are int and float, string is default.
    #Note there is also a default parameter.
    parser.add_option('-d','--dir',dest="hmm_fit_dir",type="string")
    parser.add_option('-o','--out',dest="out_path",type="string")
    parser.add_option('-t','--thresh',dest="pnathresh",type="float",default=.03)
    opts,args=parser.parse_args() #Args taken from sys.argv[1:] by default, parsed using GNU/POSIX syntax.
    if not opts.hmm_fit_dir and opts.out_path:
        parser.error("A directory for locating hmm_fit data and output file path is required.")

    print "Starting hmmprob_to_est.py with parameters:", str(opts)
    print "Free memory is %s MB" % get_free_memory()
    all_files = grab_files(opts.hmm_fit_dir)
    print "Found %s files" % len(all_files)
    d_ests = transform(all_files, opts.pnathresh)
    print "Free memory now is %s MB" % get_free_memory()
    write_csv(d_ests, opts.out_path)
    
if __name__=='__main__':
    main()