
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

import hashlib
import os
import sys

try:
    import checksums
except ImportError:
    print "Make sure you copied checksums.example.py and renamed it to checksums.py"
    sys.exit(1)

def get_dir_md5(dir_root):
    """Walk the directory and return its md5 sum
    Got some code and ideas from this post: 
    http://stackoverflow.com/questions/7325072/python-md5-of-tar-file-changes-after-file-changes-even-if-files-are-identical
    """
    hash = hashlib.md5()
    for dirpath, dirnames, filenames in os.walk(dir_root, topdown=True):
        dirnames.sort(key=os.path.normcase)
        filenames.sort(key=os.path.normcase)

        for filename in filenames:
            filepath = os.path.join(dirpath, filename)

            # If some metadata is required, add it to the checksum
            st = os.stat(filepath)

            # 1) filename
            hash.update(os.path.normcase(os.path.relpath(filepath, dir_root)))

            # 2) mtime (possibly a bad idea)
            # hash.update(struct.pack('d', st.st_mtime))
            
            # 3) size or content
            is_non_det_file = any([filename.endswith(ending) for ending in 
                checksums.non_deterministic_file_name_endings])
            if is_non_det_file:
                #just compare the file size
                #print filename,"is nondeterministic, comparing size only"
                #hash.update(bytes(st.st_size))
                pass #skip for now
            else:
                #compare content
                #print filename,"comparing content"
                f = open(filepath, 'rb')
                for chunk in iter(lambda: f.read(65536), b''):
                    hash.update(chunk)
                f.close()

    return hash.hexdigest()

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

