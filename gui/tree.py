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
Tree widget implementation.
"""
        
import Tkinter as tk
import ttk
import gui
from widgets import Scrollable

class Tree(ttk.Treeview, Scrollable):
    """ This widget acts as a wrapper around the ttk.Treeview"""
    
    def __init__(self, master, frame=True, **opts):
        """Create a new Tree, if frame=True, a container
        frame will be created and an scrollbar will be added"""
        Scrollable.__init__(self, master, ttk.Treeview, frame, **opts)
        
    def getFirst(self):
        ''' Return first selected item or None if selection empty'''
        selection = self.selection()
        if len(selection):
            return selection[0]
        return None
    
    def _moveSelection(self, moveFunc):
        item = self.selection_first()
        if item:
            item = moveFunc(item)
            if item != '':
                self.selection_set(item)
        
    def moveSelectionUp(self, e=None):
        ''' change selection to previous item '''
        self._moveSelection(self.prev)
    
    def moveSelectionDown(self, e=None):
        ''' change selection to to next item '''
        self._moveSelection(self.next)
        
    def moveItemUp(self, e=None):
        '''if selected item is not the first move up one position'''
        item = self.selection_first()
        if item:
            index = self.index(item)
            if index > 0:
                self.move(item, '', index-1)
                
    def moveItemDown(self, e=None):
        '''if selected item is not the first move up one position'''
        item = self.selection_first()
        if item:
            index = self.index(item)
            if self.next(item) != '':
                self.move(item, '', index+1)
                
    def clear(self):
        ''' remove all items '''
        childs = self.get_children('')
        for c in childs:
            self.delete(c)
            
            
class TreeProvider():
    """ Class class will serve to separete the logic of feed data
    from the graphical Tree build. Subclasses should implement 
    the abstract methods """
        
    def getColumns(self):
        """Return a list of tuples (c, w) where:
        c: is the column name and index
        w: is the column width
        """
        pass
    
    def getObjects(self):
        """Return the objects that will be inserted in the Tree"""
        pass
    
    def getObjectInfo(self, obj):
        """ This function will be called by the Tree with each
        object that will be inserted. A dictionary should be 
        returned with the possible following entries:
        'key': the key value to inserte in the Tree
        'text': text of the object to be displayed (if not passed the 'key' will be used)
        'image': image path to be displayed as icon (optional)
        'parent': the object's parent in which insert this object (optional)
        """
        pass
        
        
class BoundTree(Tree):
    """ This class is base on Tree but fetch the
    items from a TreeProvider, which provides columns
    values for each item and items info to insert into the Tree """
    def __init__(self, master, provider, frame=True, **opts):
        """Create a new Tree, if frame=True, a container
        frame will be created and an scrollbar will be added"""
        # Get colums to display and width
        cols = provider.getColumns()
        colsTuple = tuple([c[0] for c in cols[1:]])
        Tree.__init__(self, master, frame, columns=colsTuple, **opts)
        # Set the special case of first tree column
        self.heading('#0', text=cols[0][0])
        self.column('#0', width=cols[0][1])
        # Set other columns
        for c, w in cols[1:]:
            self.column(c, width=w)
            self.heading(c, text=c)
        self.grid(row=0, column=0, sticky='news')
        self._objects = provider.getObjects()
        for obj in self._objects:
            objDict = provider.getObjectInfo(obj)
            parent = objDict.get('parent', None)
            if parent is None:
                parentId = ''
            else:
                parentId = parent._treeId # Previously set
            key = objDict.get('key')
            text = objDict.get('text', key)
            image = objDict.get('image', '')
            if len(image):
                image = self.getImage(image)
            values = objDict.get('values', ())
            obj._treeId = self.insert(parentId, 'end', objDict['key'],
                        text=text, image=image, values=values)
        
        