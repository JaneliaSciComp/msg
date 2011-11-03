
"""
Python utility functions used in MSG.

"""

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

