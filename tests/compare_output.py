
"""
Run a md5sum on the directories listed in the checksums file.
Output which directories don't match their checksums.

This expects absolute paths in checksums.py.

The basic work flow for these tests is:

Do one good run with expected succesful data
if the run is succesful and the data looks good:
    Update checksums.py with the paths to the known good output folders and leave checksums as ''
    Run this script to get the good checksums.
    Update checksums.py with these.
    Edit MSG code and run it.
    After each run, run this script to check if output has deviated.

Obviously if you make changes to the code that would alter the output data then this
script isn't for you.

"""

import tarfile
import hashlib
import os
import sys

try:
    import checksums
except ImportError:
    print "Make sure you copied checksums.example.py and renamed it to checksums.py"
    sys.exit(1)

def get_dir_md5(dir_path):
    """Build a tar file of the directory and return its md5 sum"""
    temp_tar_path = 'msg_tests.tar'
    t = tarfile.TarFile(temp_tar_path,mode='w')  
    t.add(dir_path)
    t.close()

    m = hashlib.md5()
    m.update(open(temp_tar_path,'rb').read())
    ret_str = m.hexdigest()

    #delete tar file
    os.remove(temp_tar_path)
    return ret_str

def main():
    """Parse command line args, and call appropriate functions."""
    for (path, md5sum) in checksums.paths_w_checksums:
        new_md5 = get_dir_md5(path)
        if new_md5 != md5sum:
            print "MD5 sum of path '%s' has changed from '%s' to '%s'" % (path,
                md5sum, new_md5)
    print 'done'

if __name__=='__main__':
    main()

