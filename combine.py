
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
import uuid
import gc

import numpy
import numpy.lib.recfunctions

from msglib import trace, get_free_memory


# -------------- SETTINGS ----------------

# Assumes the files matches this pattern relative to hmm_fit (or other specificied directory)
GLOB_PATTERN = '/*/*-hmmprob.RData'
DEBUG = True 

# ----------------------------------------

def grab_files(dir):
    """Example from Toy data:
    glob.glob('hmm_fit/*/*-hmmprob.RData.chrom.*.csv')
    ['hmm_fit/indivF11_GTTACG/indivF11_GTTACG-hmmprob.RData.chrom.2R.csv', 'hmm_fit/indivE2_CAGCCG/indivE2_CAGCCG-hmmprob.RData.chrom.2R.csv', 
    'hmm_fit/indivG7_CTTGCG/indivG7_CTTGCG-hmmprob.RData.chrom.2R.csv', ...]
    """
    glob_pattern = dir.rstrip('/') + GLOB_PATTERN
    files = glob.glob(glob_pattern)
    print "found %s input files" % len(files)
    return files

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
    cat(sprintf('%s%%s.tsv\n', obj) )
    write.table( .x, file=paste("%s", obj, ".tsv", sep=""), sep="\t", col.names = NA,
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
    write.table( obj[[nam]], file=paste("%s", nam, ".tsv", sep=""), sep="\t", col.names = NA,
            qmethod = "double")
    cat(sprintf('%s%%s.tsv\n', nam) )
    }
   )
}
resave('%s')"""
    files_prefix = 'temp-' + str(uuid.uuid4())
    if target_object:
        r_code = highly_targeted_r_code % (target_object, files_prefix, files_prefix, rdata_file_path)
    else:
        r_code = generic_r_code % (files_prefix, files_prefix, rdata_file_path)
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
                dtype={'names': ('pos', 'par1', 'par2'), 'formats': ('a100', 'f8', 'f8')}
                )
            os.remove(csv_path)
            yield array, indiv, csv_path.replace(files_prefix,'').strip('.tsv')

def input_data_sets(dir):
    for path in grab_files(dir):
        for (array, ind, chrom) in rdata_to_numpy_arrays(path, target_object = 'dataa'):
            yield array, ind, chrom

@trace
def merge(dir):
    """
    Combine all individuals and datapoints with one row per individual, with columns
    being chrom:position.  Interpolate missing values in some cases.  (The R code
    that we're trying to replicate was funny with this so there are a few special cases,
    see code)
    Write out one tsv file for each parent.
    """
    
    #Combine all individuals/positions into a big dictionary (think of it like a sparse table)
    #for each parent
    dp1, dp2 = {}, {}
    for (array, ind, chrom) in input_data_sets(dir):
        print ind, chrom, len(array), "records"
        for x in array:
            key = (ind, chrom, int(x['pos']))
            dp1[key] = x['par1']
            dp2[key] = x['par2']

    gc.collect()
    print "Done loading rdata files."
    print "Free memory is %s MB" % get_free_memory()    

    #write out to files and interpolate as we go.  The R code we're replacing had some weird special cases so look out for those.
    for (fname, dp) in (('ancestry-probs-par1.tsv',dp1),('ancestry-probs-par2.tsv',dp2)):
        if DEBUG:
            fname = 'test.' + fname
        print "Compiling data for file",fname
        #Get all positions (chrom,pos) sorted by chrom, then by position
        positions = sorted(set([(k[1],k[2]) for k in dp.keys()]))
        header = [''] + [''.join((p[0],':',str(p[1]))) for p in positions]
        #Get all individuals, sorted
        inds = sorted(set([k[0] for k in dp.keys()]))
        #Build up each row to be written to the file (all individuals x all positions)
        outrows = []
        for ind in inds:
            print "    ",ind
            #initialize/clear out bookkeeping variables
            last_pos_w_val, last_val, last_chrom, to_interpolate = None, None, None, []

            outrow = [ind] #first column is individual name
            for (chrom,pos) in positions:
                
                # Handle switching to new chromosome
                if chrom != last_chrom:
                    #set any positions waiting for interpolation to 0 since we've reached the end of the chrom
                    #however we wan't to leave as NA and not interpolate between last_pos_w_val and end of chrom
                    #because that's what R did.
                    for (update_pos, insert_loc) in to_interpolate:
                        if update_pos < last_pos_w_val:
                            outrow[insert_loc] = "0"
                    #clear out bookkeeping vars on new chrom
                    last_pos_w_val, last_val, last_chrom, to_interpolate = None, None, None, []

                key = (ind,chrom,pos)
                if (key in dp) and ((dp[key]>.0000005) or (last_val and last_val >.0000005)):
                    # This condition is checking if A. data exists for this position and it's non-zero OR B. data exists and the last value seen was non-zero.
                    # These are cases were we want to use this value and last seen value to interpolate positions in the interpolation queue.
                    
                    # Store value in outrow to be written to file
                    outrow.append("%.6f" % round(dp[key],6))
                    #interpolate any positions waiting for a new value
                    for (update_pos, insert_loc) in to_interpolate:
                        if update_pos < last_pos_w_val:
                            outrow[insert_loc] = "0" # zero out any pending positions before the last value we saw since this is what R did.
                        else:
                            insert_val = last_val + ((dp[key] - last_val) * (float(update_pos - last_pos_w_val) / (pos - last_pos_w_val)))
                            outrow[insert_loc] = "%.6f" % round(insert_val,6)
                    to_interpolate = [] #since all pending positions have been interpolated, clear this out
                    last_pos_w_val, last_val = pos, dp[key]
                elif last_val and not (key in dp):
                    #If a value has been seen for this chrom, we'll want to start interpolating
                    #Add a placeholder to outrow
                    outrow.append('NA') #
                    #Mark position for later interpolation
                    to_interpolate.append((pos, len(outrow) - 1))
                else:
                    #don't interpolate
                    if key in dp:
                        #data exists for key but it's 0, Store value in outrow, but update bookkeeping vars
                        outrow.append("%.6f" % round(dp[key],6)) #should be 0
                        #still count 0 as a last value for interpolation
                        last_pos_w_val, last_val = pos, dp[key]
                    else:
                        outrow.append('NA')
                last_chrom = chrom

            #set any positions waiting for interpolation to 0 since we've reached the end of the individual
            #however we wan't to leave as NA and not interpolate between last_pos_w_val and end
            #because that's what R did.
            for (update_pos, insert_loc) in to_interpolate:
                if update_pos < last_pos_w_val:
                    outrow[insert_loc] = "0"

            outrows.append(outrow)
        print "Writing file",fname
        csvout = csv.writer(open(fname,'wb'), delimiter='\t', quoting=csv.QUOTE_MINIMAL)
        csvout.writerow(header)
        csvout.writerows(outrows)
        gc.collect()

@trace
def main():
    """Parse command line args, and call appropriate functions."""
    #disable garbage collection for a 10% speed boost
    gc.disable()
    
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
    