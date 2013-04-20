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
        
import Tkinter as tk

import pyworkflow.gui as gui
from pyworkflow.object import *
from pyworkflow.mapper import SqliteMapper, XmlMapper
from pyworkflow.gui.tree import BoundTree, DbTreeProvider
from pyworkflow.protocol import *
from pyworkflow.protocol.params import *
from pyworkflow.tests.tester import *
from pyworkflow.em import *

            
def getMapper(fn, classesDict):
    """Select what Mapper to use depending on
    the filename extension"""
    if fn.endswith('.xml'):
        return XmlMapper(fn, classesDict)
    elif fn.endswith('.sqlite'):
        return SqliteMapper(fn, classesDict)
    return None

    
class BrowserWindow(gui.Window):
    def __init__(self, title, provider, master=None, **args):
        if 'minsize' not in args:
            args['minsize'] = (800, 400)
        gui.Window.__init__(self, title, master, **args)

        tree = BoundTree(self.root, provider)
        #self.populateTree(tree, objList)
        tree.grid(row=0, column=0, sticky='news')
        self.itemConfig = tree.itemConfig

    def populateWithObject(self, tree, prefix, obj):
        
        sid = obj.getName()
        if len(sid) == 0:
            sid = obj.strId()
        item = tree.insert(prefix, 'end', sid, text=t, values=(cls,))
        haschilds = False

        if not haschilds:
            img = self.getImage('step.gif')
            tree.item(item, image=img)
            
        return item
            
    def populateTree(self, tree, objList):
        """Populate a tree from an object list"""
        for obj in objList:
            item = self.populateWithObject(tree, '', obj)
            tree.item(item, open=True)
    
if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        path = sys.argv[1]
        provider = DbTreeProvider(path, globals())
        window = BrowserWindow("Browsing: " + path, provider)    
        window.show()
    else:
        print "usage: pw_browser.py SQLITE_DB"
