
"""
Combines each individual's hmmprob.RData file into two summary files
(linearly interpolating missing values)


Usage:
python msg/combine.py
msg/combine.py -d /groups/stern/home/sternd/svb_mausec2/hmm_fit
"""

import os
import sys
import csv
import glob
import optparse
import subprocess

import numpy
import numpy.lib.recfunctions
#import rpy2.robjects as robjects
#robjects.r['load'](".RData")
#robjects.r['y']

from msglib import trace, get_free_memory


# -------------- SETTINGS ----------------

# Assumes the files matches this pattern relative to hmm_fit (or other specificied directory)
GLOB_PATTERN = '/*/*-hmmprob.RData'


# ----------------------------------------

def grab_files(dir):
    """Example from Toy data:
    glob.glob('hmm_fit/*/*-hmmprob.RData.chrom.*.csv')
    ['hmm_fit/indivF11_GTTACG/indivF11_GTTACG-hmmprob.RData.chrom.2R.csv', 'hmm_fit/indivE2_CAGCCG/indivE2_CAGCCG-hmmprob.RData.chrom.2R.csv', 
    'hmm_fit/indivG7_CTTGCG/indivG7_CTTGCG-hmmprob.RData.chrom.2R.csv', ...]
    """
    glob_pattern = dir.rstrip('/') + GLOB_PATTERN
    return glob.glob(glob_pattern)

def parse_path(path):
    """Get ind name, and chrom from file path"""
    dir, filename = os.path.split(path)
    name_parts = filename.split('.')
    ind_name = filename.split('-hmmprob')[0]
    return ind_name

def rdata_to_numpy_arrays(rdata_file_path, target_object=None):
    """Call out to R (on $PATH) to convert rdata file to one or more
    CSV files.  Load CSV files into numpy arrays and delete CSV files.
    
    If target_object is None, it will try to use generic code to find all the dataframes.  Otherwise
    it will try to home in on the target_object and find the frames within that.
    
    (R code from http://stackoverflow.com/questions/13189467/how-to-convert-rdata-format-into-text-file-format)
    More discussion here:
    http://stackoverflow.com/questions/23413728/converting-rdata-files-to-csv-error-in-data-frame-arguments-imply-differing-nu
    """    
    generic_r_code = """\
resave <- function(file){
  e <- new.env(parent = emptyenv())
  load(file, envir = e)
  objs <- ls(envir = e, all.names = TRUE)
  for(obj in objs) {
    .x <- get(obj, envir =e)
    cat(sprintf('%%s.tsv\n', obj) )
    write.table( .x, file=paste(obj, ".tsv", sep=""), sep="\t", col.names = NA,
            qmethod = "double")    
  }
}
resave('%s')"""

    highly_targeted_r_code = """\
resave <- function(file){
  e <- new.env(parent = emptyenv())
  load(file, envir = e)
  obj <- get('%s', envir =e)
  lapply( names(obj), function(nam) {
    write.table( obj[[nam]], file=paste(nam, ".tsv", sep=""), sep="\t", col.names = NA,
            qmethod = "double")
    cat(sprintf('%%s.tsv\n', nam) )
    }
   )
}
resave('%s')"""
    if target_object:
        r_code = highly_targeted_r_code % (target_object, rdata_file_path)
    else:
        r_code = generic_r_code % (rdata_file_path)
    #print r_code
    command = ["Rscript","-","-"] #"-" to tell Rscript to use pipes for input and output
    #print ' '.join(command)
    rscript = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    file_list = rscript.communicate(r_code)[0]
    indiv = parse_path(rdata_file_path)
    for csv_path in file_list.splitlines():
        if csv_path.lower().endswith('.tsv'):
            #Note: Setting the comments parameter below is a numpy hack to make it not look
            #for comments in our data file
            array = numpy.loadtxt(csv_path, skiprows=1, usecols=(1,21,23), delimiter="\t",
                comments="wewillneverseethisstringinafile15",
                dtype={'names': ('pos', indiv+'-par1', indiv+'-par2'), 'formats': ('a100', 'f8', 'f8')})
            os.remove(csv_path)
            yield array, indiv, csv_path.strip('.tsv')

def input_data_sets(dir):
    for path in grab_files(dir):
        for (array, ind, chrom) in rdata_to_numpy_arrays(path, target_object = 'dataa'):
            yield array, ind, chrom

@trace
def merge(dir):
    """
    This is my attempt at using numpy for the join, Ultimately it seemed to insist on
    creating new columns even if the field names were the same.
    Here's some discussion on the issues:
    http://stackoverflow.com/questions/23500754/numpy-how-to-outer-join-arrays

    array will look like this:
    
array([('2:6506', 4.6725971801473496e-25, 0.99999999995088695),
       ('2:6601', 2.2452745388799898e-27, 0.99999999995270605),
       ('2:21801', 1.9849650921836601e-31, 0.99999999997999001), ...,
       ('2:45164194', 1.0413482803123399e-24, 0.99999999997453404),
       ('2:45164198', 1.09470356446595e-24, 0.99999999997635303),
       ('2:45164519', 3.7521365799080699e-24, 0.99999999997453404)], 
      dtype=[('pos', '|S100'), ('indiv8.3.1_4_tatccTCAAGCCTCCTCGCACtga-par1', '<f8'), ('indiv8.3.1_4_tatccTCAAGCCTCCTCGCACtga-par2', '<f8')])

    """
    final_array = None
    
    for (array, ind, chrom) in input_data_sets(dir):
        #Put chrom in pos
        array['pos'] = [''.join((chrom,':',x)) for x in array['pos']]
        print ind, chrom, len(array), "records"
        
        if final_array is None:
            final_array = array
            continue
            
        final_array = numpy.lib.recfunctions.join_by('pos', final_array, array, jointype='outer',
            r1postfix='', r2postfix='')
    final_array.fill_value = numpy.NaN
    final_array = final_array.filled()

    #TODO: separate out par1, par2, transpose final_array and write out to file


@trace
def main():
    """Parse command line args, and call appropriate functions."""
    usage="""\nusage: %prog [options]\n"""
    parser = optparse.OptionParser(usage=usage)
    #Other option types are int and float, string is default.
    #Note there is also a default parameter.
    parser.add_option('-d','--dir',dest="hmm_fit_dir",type="string")
    #?? Need these ??  -c $params{'chroms'} -p $params{'chroms2plot'} -d hmm_fit -t $params{'thinfac'} -f $params{'difffac'} -b $params{'barcodes'} -n $params{'pnathresh'}
    #parser.add_option('-o','--out',dest="out_path",type="string")
    #parser.add_option('-t','--thresh',dest="pnathresh",type="float",default=.03)
    opts,args=parser.parse_args() #Args taken from sys.argv[1:] by default, parsed using GNU/POSIX syntax.
    if not opts.hmm_fit_dir:
        parser.error("A directory for locating hmm_fit data is required.")

    print "Starting combine.py with parameters:", str(opts)
    print "Free memory is %s MB" % get_free_memory()
    merge(opts.hmm_fit_dir)
    
if __name__=='__main__':
    main()