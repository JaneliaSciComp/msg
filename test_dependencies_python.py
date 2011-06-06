#!/usr/env python
import sys


try:
     import pysam 
     import Bio
except ImportError, e:
    print >>sys.stderr, "Required module not found: %s" % e
    sys.exit(2)

print("All required python modules found")
sys.exit(0)
