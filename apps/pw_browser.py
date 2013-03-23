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
from gui.gui import getImage
"""
Object browser
"""
        
import Tkinter as tk
from pyworkflow.object import *

from pyworkflow.mapper import SqliteMapper, XmlMapper
import gui
from gui.widgets import Tree
from protocol import *
from protocol.params import *

class MyStep(Step):
    def __init__(self):
        Step.__init__(self)
        self.defineInputs(x=Integer(1), y=Float(2), z=String("abc"), b=Boolean(True))
        
    def __str__(self):
        s = ''
        for k, v in self.getAttributesToStore():
            s += '%s = %s\n' % (k, str(v))
        return s

    def hasValue(self):
        return True  

class Complex(Object):
    def __init__(self, imag=0., real=0., **args):
        Object.__init__(self, **args)
        self.imag = Float(imag)
        self.real = Float(real)
        
    def __str__(self):
        return '(%s, %s)' % (self.imag, self.real)
    
    def __eq__(self, other):
        return self.imag == other.imag and \
            self.real == other.real
            
    def hasValue(self):
        return True
    

def populateWithObject(tree, prefix, obj):
    cls = obj.getClassName()
    if obj.parent_id is None:
        t = cls
    else:
        t = obj.name.split('.')[-1] 
        if  t.startswith('__item__'):
            t = "%s [%s]" % (cls, t.replace('__item__', ''))
    if obj.value:
        t += " = %s" % str(obj.value)
    item = tree.insert(prefix, 'end', obj.name, text=t, values=(cls,))
    haschilds = False
    for k, v in obj.getAttributesToStore():
        populateWithObject(tree, obj.name, v)
        haschilds = True
    if not haschilds:
        img = getImage('step.gif')
        tree.item(item, image=img)
        
    return item
        
def populateTree(tree, objList):
    """Populate a tree from an object list"""
    for obj in objList:
        item = populateWithObject(tree, '', obj)
        tree.item(item, open=True)
    
def getMapper(fn, classesDict):
    """Select what Mapper to use depending on
    the filename extension"""
    if fn.endswith('.xml'):
        return XmlMapper(fn, classesDict)
    elif fn.endswith('.sqlite'):
        return SqliteMapper(fn, classesDict)
    return None

if __name__ == '__main__':
    import sys
    fn = sys.argv[1]
    mapper = getMapper(fn, globals())
    objList = mapper.getAll()
    
    window = gui.Window("Object Browser")
    
    columns = ('Class', )
    c = columns[0]
    tree = Tree(window.root, columns=columns)
    tree.heading(c, text=c) 
    tree.heading('#0', text='Object')
    tree.column('#0', minwidth=300)
    populateTree(tree, objList)
    tree.grid(row=0, column=0, sticky='news')
    
    window.show()
