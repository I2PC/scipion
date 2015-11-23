#!/usr/bin/env python


import os
import sys

from pyworkflow.em.viewer import DataView
from pyworkflow.gui.tree import DbTreeProvider
from pyworkflow.gui.browser import FileBrowserWindow
            
   
    
def showDir(path):
    window = FileBrowserWindow("Browsing: " + path, path=path) 
    window.show()
    
def showFile(path, viewParams):
    DataView(path, viewParams).show()
        


if __name__ == '__main__':    
    
    args = {'-i': []}
    lastOpt = '-i'
    
    if len(sys.argv) > 1 or '-h' in sys.argv or '--help' in sys.argv:
        for a in sys.argv[1:]:
            if a.startswith('-'):
                lastOpt = a
                if lastOpt not in args:
                    args[lastOpt] = []
            else:
                args[lastOpt].append(a)
                
        inputFiles = args['-i']
        del args['-i']
        viewParams = {}
        for k, v in args.iteritems():
            viewParams[k.replace('-', '')] = ' '.join(v)
        
        for fn in inputFiles:
            if os.path.isdir(fn):
                showDir(fn)
            else:
                showFile(fn, viewParams)
    else:
        print "usage: pw_viewer.py file1 [file2 file3 ... fileN]"
     
