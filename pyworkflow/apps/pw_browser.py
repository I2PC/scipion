#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
Object browser
"""
        
import sys
import Tkinter as tk

from pyworkflow.gui.tree import DbTreeProvider
from pyworkflow.gui.browser import *
            
            
class EMTreeProvider(DbTreeProvider):
    """Retrieve the elements from the database"""
    def __init__(self, dbName):
        print 'complex in globals: ', 'Complex' in globals()
        DbTreeProvider.__init__(self, dbName, globals())
        self.viewer = XmippViewer()
        
    def show(self, obj):
        self.viewer.visualize(obj)
        
    def getObjectPreview(self, obj):
        desc = "<name>: " + obj.getName() + "\n " + str(obj)
        
        return (None, desc)
    
    
        
USAGE = "usage: pw_browser.py [db|dir] path"

if __name__ == '__main__':

    if len(sys.argv) == 3:
        browseMode = sys.argv[1]
        path = sys.argv[2]
        if browseMode == 'dir':
            window = FileBrowserWindow("Browsing: " + path, path=path) 
        elif browseMode == 'db':
            window = BrowserWindow("Browsing: " + path)    
            provider = EMTreeProvider(path)
            browser = ObjectBrowser(window.root, provider)
            window.setBrowser(browser)
        else:
            print "Unknown mode %s\n%s" % (browseMode, USAGE)
        window.show()
    else:
        print USAGE

