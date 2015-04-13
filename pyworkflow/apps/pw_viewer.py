#!/usr/bin/env python


import sys
from pyworkflow.em.viewer import DataView



if __name__ == '__main__':    
    
    if len(sys.argv) > 1:
        for fn in sys.argv[1:]:
            DataView(fn).show()
    else:
        print "usage: pw_viewer.py file1 [file2 file3 ... fileN]"
     
