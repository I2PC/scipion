#!/usr/bin/env python


import os
import sys

from pyworkflow.em.viewer import DataView
from pyworkflow.gui.tree import DbTreeProvider
from pyworkflow.gui.browser import FileBrowserWindow
            
   
    
def showDir(path):
    window = FileBrowserWindow("Browsing: " + path, path=path) 
    window.show()
    
    
def showSqlite(path):
    # FIXME: Develop a way to display sqlite files
    window = BrowserWindow("Browsing: " + path)    
    provider = EMTreeProvider(path)
    browser = ObjectBrowser(window.root, provider)
    window.setBrowser(browser)
    window.show()
        

def showFile(path):
    DataView(path).show()
        


if __name__ == '__main__':    
    
    if len(sys.argv) > 1:
        for fn in sys.argv[1:]:
            if os.path.isdir(fn):
                showDir(fn)
            else:
                showFile(fn)
    else:
        print "usage: pw_viewer.py file1 [file2 file3 ... fileN]"
     
