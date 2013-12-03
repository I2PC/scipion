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

from pyworkflow.object import Scalar
from pyworkflow.mapper import SqliteMapper
import gui
from widgets import Scrollable

class Tree(ttk.Treeview, Scrollable):
    """ This widget acts as a wrapper around the ttk.Treeview"""
    _images = {}
    
    def __init__(self, master, frame=True, **opts):
        """Create a new Tree, if frame=True, a container
        frame will be created and an scrollbar will be added"""
        Scrollable.__init__(self, master, ttk.Treeview, frame, **opts)
        
    def getImage(self, img):
        return gui.getImage(img, Tree._images)
    
    def getFirst(self):
        """ Return first selected item or None if selection empty"""
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
        """ change selection to previous item """
        self._moveSelection(self.prev)
    
    def moveSelectionDown(self, e=None):
        """ change selection to to next item """
        self._moveSelection(self.next)
        
    def moveItemUp(self, e=None):
        """if selected item is not the first move up one position"""
        item = self.selection_first()
        if item:
            index = self.index(item)
            if index > 0:
                self.move(item, '', index-1)
                
    def moveItemDown(self, e=None):
        """if selected item is not the first move up one position"""
        item = self.selection_first()
        if item:
            index = self.index(item)
            if self.next(item) != '':
                self.move(item, '', index+1)
                
    def clear(self):
        """ remove all items """
        childs = self.get_children('')
        for c in childs:
            self.delete(c)
            
    def selectItem(self, index):
        """ Select the item at the position index """
        child = self.get_children('')[index]
        self.selection_set(child)
            
            
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
        'key': the key value to insert in the Tree
        'text': text of the object to be displayed (if not passed the 'key' will be used)
        'image': image path to be displayed as icon (optional)
        'parent': the object's parent in which insert this object (optional)
        """
        pass
    
    def getObjectPreview(self, obj):
        """ Should return a tuple (img, desc),
        where img is the preview image and 
        desc the description string. 
        """
        return (None, None)
    
    def getObjectActions(self, obj):
        """ Return a list of tuples (key, action)
        were keys are the string
        options that will be display in the context menu
        and the actions are the functions to call when
        the specific action is selected.
        The first action in the list will be taken
        as the default one when the element is double-clicked.
        """
        return []
        
        
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
        
        self.menu = tk.Menu(self, tearoff=0)
        self.setProvider(provider)
        
        self.bind("<Button-3>", self._onRightClick)
        # Hide the right-click menu
        self.bind('<FocusOut>', self._unpostMenu)
        self.bind("<Key>", self._unpostMenu)
        self.bind('<Button-1>', self._unpostMenu)
        self.bind('<Double-1>', self._onDoubleClick)
        self.bind('<<TreeviewSelect>>', self._onClick)
        
    def setProvider(self, provider):
        """ Set new provider and updated items. """
        self.provider = provider
        self.update()
        
        
    def _unpostMenu(self, e=None):
        self.menu.unpost()
        
    def _onClick(self, e=None):
        if hasattr(self, 'itemClick'):
            obj = self._objDict[self.getFirst()]
            self.itemClick(obj)
            
    def _onDoubleClick(self, e=None):  
        obj = self._objDict[self.getFirst()]
        if hasattr(self, 'itemDoubleClick'):
            self.itemDoubleClick(obj)
        else: # If not callback, use default action
            actions = self.provider.getObjectActions(obj)
            if len(actions):
                actions[0][1]() # actions[0] = first Action, [1] = the action callback
            
    def _onRightClick(self, e=None):
        item = self.identify('item', e.x, e.y)
        unpost = True
        if len(item):
            self.selection_set(item)
            obj = self._objDict[item]
            actions = self.provider.getObjectActions(obj)
            if len(actions):
                self.menu.delete(0, tk.END)
                for a in actions:
                    if a is None: 
                        self.menu.add_separator()
                    else:
                        img = ''
                        if len(a) > 2: # image for the action
                            img = self.getImage(a[2])
                        self.menu.add_command(label=a[0], command=a[1], 
                                              image=img, compound=tk.LEFT)
                self.menu.post(e.x_root, e.y_root)
                unpost = False
        if unpost:
            self._unpostMenu()        
                    
    def update(self):
        self.clear()
        self._objDict = {} # Store the mapping between Tree ids and objects
        self._objects = self.provider.getObjects()
        for obj in self._objects:
            # If the object is a pointer that has a null value do not show
            #if ((not obj.isPointer()) or (obj.isPointer() and obj.get() is not None)): 
            objDict = self.provider.getObjectInfo(obj)
            if objDict is not None:
                key = objDict.get('key')
                text = objDict.get('text', key)
                parent = objDict.get('parent', None)
                open = objDict.get('open', False)
                
                if parent is None:
                    parentId = ''
                else:
                    if hasattr(parent, '_treeId'): # This should happens always
                        parentId = parent._treeId # Previously set
                    else:
                        parentId = ''
                        text += '---> Error: parent not Inserted'
                image = objDict.get('image', '')
                if len(image):
                    image = self.getImage(image)
                    if image is None:
                        image = ''
                values = objDict.get('values', ())
                obj._treeId = self.insert(parentId, 'end', key,
                            text=text, image=image, values=values)
                self._objDict[obj._treeId] = obj
                if open:
                    self.itemConfig(obj, open=open)

    def itemConfig(self, obj, **args):
        """ Configure inserted items. """
        self.item(obj._treeId, **args)
        
              
class ObjectTreeProvider(TreeProvider):
    """ Populate Tree from Objects. """
    def __init__(self, objList=None):
        self.objList = objList
        self.getColumns = lambda: [('Object', 300), ('Id', 70), ('Class', 150)]
        self._parentDict = {}
    
    def getObjectInfo(self, obj):
        if obj.isPointer() and not obj.hasValue():
            return None
        cls = obj.getClassName()
        if obj.getName() is None:
            t = cls
        else:
            t = obj.getName().split('.')[-1] 
            if  t.startswith('__item__'):
                t = "%s [%s]" % (cls, t.replace('__item__', ''))
        if obj.get() is not None:
            t += " = %s" % str(obj)
            
        info = {'key': obj.getObjId(), 'parent': self._parentDict.get(obj.getObjId(), None),
                'text': t, 'values': (obj.strId(), cls)}
        if issubclass(obj.__class__, Scalar):
            info['image'] = 'step.gif'
            
        return info
    
    def getObjectPreview(self, obj):
        return (None, None)
    
    def getObjectActions(self, obj):
        return []
        
    def _getObjectList(self):
        """Retrieve the object list"""
        return self.objList
    
    def getObjects(self):
        objList = self._getObjectList()
        self._parentDict = {}
        childs = []
        for obj in objList:
            childs += self._getChilds(obj)
        objList += childs
        return objList
        
    def _getChilds(self, obj):
        childs = []
        grandchilds = []
        for _, v in obj.getAttributesToStore():
            childs.append(v)
            self._parentDict[v.getObjId()] = obj
            grandchilds += self._getChilds(v)
        childs += grandchilds
        return childs

    
class DbTreeProvider(ObjectTreeProvider):
    """Retrieve the elements from the database"""
    def __init__(self, dbName, classesDict):
        ObjectTreeProvider.__init__(self)
        self.mapper = SqliteMapper(dbName, classesDict)
    
    def _getObjectList(self):
        return self.mapper.selectAll()
    
    
class ProjectRunsTreeProvider(TreeProvider):
    """ Provide run list from a project
    to populate a tree.
    """
    def __init__(self, project):
        self.project = project
        self._objDict = {}
    
    def getObjects(self):
        return self.project.getRuns() 
        
    def getColumns(self):
        return [('Run', 250), ('State', 100), ('Time', 100)]
    
    def getObjectInfo(self, obj):
        objId = obj.getObjId()
        self._objDict[objId] = obj
        
        info = {'key': obj.getObjId(),
                'text': obj.getRunName(),
                'values': (obj.getStatusMessage(), obj.getElapsedTime())}
        objPid = obj.getObjParentId()
        if objPid in self._objDict:
            info['parent'] = self._objDict[objPid]
      
        return info

