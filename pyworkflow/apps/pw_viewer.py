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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

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
    
    if '-h' in sys.argv or '--help' in sys.argv:
        print "usage: scipion view [file1 file2 file3 ... fileN]"
        
    else:
        if len(sys.argv) == 1: # no extra arguments, show current directory
            showDir(os.getcwd())
        else:
            args = {'-i': []}
            lastOpt = '-i'
            
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
     
