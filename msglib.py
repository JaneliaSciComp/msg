
"""
Python utility functions used in MSG.

"""

from Queue import Queue
from threading import Thread

class bcolors(object):
   OKMAGENTA = '\033[95m'
   OKBLUE = '\033[94m'
   OKGREEN = '\033[92m'
   WARN = '\033[93m'
   FAIL = '\033[91m'
   ENDC = '\033[0m'

   def disable(self):
       self.HEADER = ''
       self.OKBLUE = ''
       self.OKGREEN = ''
       self.WARNING = ''
       self.FAIL = ''
       self.ENDC = ''

def trace(fn):
    """A decorator to time your functions"""
    from time import time
    import sys
    def trace_func(*args, **kwargs):
        print fn.__name__ + '...',
        sys.stdout.flush()
        beg = time()
        ret = fn(*args, **kwargs)
        tot = time() - beg
        print '%.3f' % tot
        return ret
    return trace_func

def sort_unique(file_path):
    """
    open a file, remove duplicates (by first column) keeping last occurance
    and sort by first column.

    Note: I also tried calling sort -u as an alternative [1] to this function
    however when that dedupes it keeps the first occurance of a duplicate
    and we want the last.  It also seemed slower.
    [1]
    #sorts only on first column, dedupes only by first column
    out = commands.getoutput('sort -ug -o "%s" "%s"' % (file_path, file_path))
    """
    lines = open(file_path,'r').readlines()
    outfile = open(file_path,'w')
    rows = [line.strip().split('\t') for line in lines]
    drows = {}
    for row in rows:
        drows[int(row[0])]='\t'.join(row[1:])
    rows = sorted(drows.items())
    kv_items = ["\t".join(map(str, kv)) for kv in rows]
    outfile.write('\n'.join(kv_items))
    outfile.write('\n')                    
    outfile.close()

