# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
This modules implements the automatic
creation of protocol form GUI from its
params definition.
"""
import os
import Tkinter as tk
import ttk

from pyworkflow.utils.path import getHomePath
from pyworkflow.utils.properties import Message, Icon, Color
from pyworkflow.viewer import DESKTOP_TKINTER
import gui
from gui import configureWeigths, Window
from browser import FileBrowserWindow
from widgets import Button, HotButton, IconButton
from pyworkflow.protocol.params import *
from dialog import showInfo, EditObjectDialog, ListDialog, askYesNo
from canvas import Canvas
from tree import TreeProvider, BoundTree
#from pyworkflow.em import findViewers

#-------------------- Variables wrappers around more complex objects -----------------------------

class BoolVar():
    """Wrapper around tk.IntVar"""
    def __init__(self, value=False):
        self.tkVar = tk.IntVar()
        self.set(value)
        self.trace = self.tkVar.trace
        
    def set(self, value):
        if value:
            self.tkVar.set(1)
        else:
            self.tkVar.set(0)    
            
    def get(self):
        return self.tkVar.get() == 1   
    
    
class PointerVar():
    """Wrapper around tk.StringVar to hold object pointers"""
    def __init__(self, protocol):
        self.tkVar = tk.StringVar()
        self.value = None
        self.trace = self.tkVar.trace
        self._protocol = protocol
        
    def set(self, value):
        self.value = value
        self.tkVar.set(self.getObjectLabel())   
     
    def getObjectLabel(self):
        """ We will try to show in the list the string representation
        that is more readable for the user to pick the desired object.
        """
        label = ''
        obj = self.value
        if obj:
            label = obj.getObjLabel()
            if not len(label.strip()):
                parent = self._protocol.mapper.getParent(obj)
                label = "%s -> %s" % (parent.getObjLabel(), obj.getLastName())
        return label
           
    def get(self):
        return self.value
    
    
class MultiPointerVar():
    """Wrapper around tk.StringVar to hold object pointers"""
    def __init__(self, provider, tree):
        self.provider = provider # keep a reference to tree provider to add or remove objects
        self.tree = tree
        self.tkVar = tk.StringVar()
        self.trace = self.tkVar.trace
        
    def set(self, value):
        if isinstance(value, Object):
            self.provider.addObject(value)
            self.tree.update()
          
    def remove(self):
        """ Remove first element selected. """
        value = self.tree.getFirst()
        if value:
            self.provider.removeObject(value)
            self.tree.update()   
        
    def get(self):
        return self.provider._objectList
    
   
class ComboVar():
    """ Create a variable that display strings (for combobox)
    but the values are integers (for the underlying EnumParam).
    """
    def __init__(self, enum):
        self.tkVar = tk.StringVar()
        self.enum = enum
        self.value = None
        self.trace = self.tkVar.trace
        
    def set(self, value):
        self.value = value
        if isinstance(value, int):
            self.tkVar.set(self.enum.choices[value])
        else:
            self.tkVar.set(value) # also support string values
                    
    def get(self):
        v = self.tkVar.get()
        self.value = None
        for i, c in enumerate(self.enum.choices):
            if c == v:
                self.value = i
            
        return self.value         
        

#---------------- Some used providers for the TREES -------------------------------

class ProtocolClassTreeProvider(TreeProvider):
    """Will implement the methods to provide the object info
    of subclasses objects(of className) found by mapper"""
    def __init__(self, protocolClassName):
        self.protocolClassName = protocolClassName
     
    def getObjects(self):
        from pyworkflow.em import findSubClasses, emProtocolsDict
        return [String(s) for s in findSubClasses(emProtocolsDict, self.protocolClassName).keys()]
        
    def getColumns(self):
        return [('Protocol', 250)]
    
    def getObjectInfo(self, obj):
        return {'key': obj.get(),
                'values': (obj.get(),)}
    

class SubclassesTreeProvider(TreeProvider):
    """Will implement the methods to provide the object info
    of subclasses objects(of className) found by mapper"""
    def __init__(self, protocol, pointerParam, selected=None):
        self.className = pointerParam.pointerClass.get()
        self.condition = pointerParam.pointerCondition.get()
        self.selected = selected
        self.protocol = protocol
        self.mapper = protocol.mapper
        
    def getObjects(self):
        return list(self.protocol.getProject().iterSubclasses(self.className, self.objFilter))
            
    def objFilter(self, obj):
        result = True
        # Do not allow to select objects that are childs of the protocol
        if (self.protocol.getObjId() and 
            self.protocol.getObjId() == obj.getObjParentId()):
            result = False
        # Check that the condition is met
        elif self.condition:
            result = obj.evalCondition(self.condition)
        return result
        
    def getColumns(self):
        return [('Object', 400), ('Info', 250)]
    
    def isSelected(self, obj):
        """ Check if an object is selected or not. """
        if self.selected:
            for s in self.selected:
                if s and s.getObjId() == obj.getObjId():
                    return True
        return False
    
    def getObjectLabel(self, obj):
        """ We will try to show in the list the string representation
        that is more readable for the user to pick the desired object.
        """
        label = obj.getObjLabel()
        if not len(label.strip()):
            parent = self.mapper.getParent(obj)
            label = "%s -> %s" % (parent.getObjLabel(), obj.getLastName())
        return label
        
    def getObjectInfo(self, obj):
        return {'key': obj.strId(), 'text': self.getObjectLabel(obj),
                'values': (str(obj),), 'selected': self.isSelected(obj)}

    def getObjectActions(self, obj):
        if isinstance(obj, Pointer):
            obj = obj.getName()
        actions = []    
        from pyworkflow.em import findViewers
        viewers = findViewers(obj.getClassName(), DESKTOP_TKINTER)
        for v in viewers:
            actions.append(('Open with %s' % v.__name__, lambda : v(project=self.protocol.getProject()).visualize(obj)))
            
        return actions
    
    
#TODO: check if need to inherit from SubclassesTreeProvider
class RelationsTreeProvider(SubclassesTreeProvider):
    """Will implement the methods to provide the object info
    of subclasses objects(of className) found by mapper"""
    def __init__(self, protocol, relationParam, selected=None):
        self.mapper = protocol.mapper
        self.selected = selected
        item = protocol.getAttributeValue(relationParam.getAttributeName())
        direction = relationParam.getDirection()
        if item is not None:
            self.getObjects = lambda: protocol.getProject().getRelatedObjects(relationParam.getName(), item, direction)
        else:
            self.getObjects = lambda: [] 
            

class MultiPointerTreeProvider(TreeProvider):
    """Store several objects to be used for Multipointer"""
    def __init__(self):
        self._objectList = []
        
    def addObject(self, obj):
        if not obj in self._objectList:
            self._objectList.append(obj)
        
    def removeObject(self, objId):
        objId = int(objId)
        for obj in self._objectList:
            if objId == obj.getObjId():
                self._objectList.remove(obj)
                break 
     
    def getObjects(self):
        return self._objectList
        
    def getColumns(self):
        return [('Object', 250)]
    
    def getObjectInfo(self, obj):
        return {'key': obj.getObjId(), 'text': (obj.getNameId(),)}  
            
   
#---------------------- Other widgets ----------------------------------------
# http://tkinter.unpythonic.net/wiki/VerticalScrolledFrame

class VerticalScrolledFrame(tk.Frame):
    """A pure Tkinter scrollable frame that actually works!
    * Use the 'interior' attribute to place widgets inside the scrollable frame
    * Construct and pack/place/grid normally
    * This frame only allows vertical scrolling

    """
    def __init__(self, parent, *args, **kw):
        tk.Frame.__init__(self, parent, *args, **kw)            

        # create a canvas object and a vertical scrollbar for scrolling it
        vscrollbar = tk.Scrollbar(self, orient=tk.VERTICAL)
        vscrollbar.pack(fill=tk.Y, side=tk.RIGHT, expand=tk.FALSE)
        canvas = Canvas(self, bd=0, highlightthickness=0,
                        yscrollcommand=vscrollbar.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=tk.TRUE)
        vscrollbar.config(command=canvas.yview)

        # reset the view
        canvas.xview_moveto(0)
        canvas.yview_moveto(0)

        # create a frame inside the canvas which will be scrolled with it
        self.interior = interior = tk.Frame(canvas)
        interior_id = canvas.create_window(0, 0, window=interior,
                                           anchor=tk.NW)

        # track changes to the canvas and frame width and sync them,
        # also updating the scrollbar
        def _configure_interior(event):
            # update the scrollbars to match the size of the inner frame
            size = (interior.winfo_reqwidth(), interior.winfo_reqheight())
            canvas.config(scrollregion="0 0 %s %s" % size)
            if interior.winfo_reqwidth() != canvas.winfo_width():
                # update the canvas's width to fit the inner frame
                canvas.config(width=interior.winfo_reqwidth())
        interior.bind('<Configure>', _configure_interior)

        def _configure_canvas(event):
            if interior.winfo_reqwidth() != canvas.winfo_width():
                # update the inner frame's width to fill the canvas
                canvas.itemconfigure(interior_id, width=canvas.winfo_width())
        canvas.bind('<Configure>', _configure_canvas)


class SectionFrame(tk.Frame):
    """This class will be used to create a frame for the Section
    That will have a header with red color and a content frame
    with white background
    """
    def __init__(self, master, label, callback=None, height=15, **args):
        headerBgColor = args.get('headerBgColor', gui.cfgButtonBgColor)
        if 'headerBgColor' in args:
            del args['headerBgColor']
        self.height = height
        tk.Frame.__init__(self, master, bg='white', **args)
        self._createHeader(label, headerBgColor)
        self._createContent()
        
    def _createHeader(self, label, bgColor):
        self.headerFrame = tk.Frame(self, bd=2, relief=tk.RAISED, bg=bgColor)
        self.headerFrame.grid(row=0, column=0, sticky='new')
        configureWeigths(self.headerFrame)
        self.headerFrame.columnconfigure(1, weight=1)
        #self.headerFrame.columnconfigure(2, weight=1)
        self.headerLabel = tk.Label(self.headerFrame, text=label, fg='white', bg=bgColor)
        self.headerLabel.grid(row=0, column=0, sticky='nw')
        
    def _createContent(self):
        canvasFrame = tk.Frame(self, bg='white')
        configureWeigths(self, row=1)
        configureWeigths(canvasFrame)
        self.canvas = Canvas(canvasFrame, width=625, height=self.height, bg='white') 
        self.canvas.grid(row=0, column=0, sticky='news')
        canvasFrame.grid(row=1, column=0, sticky='news')
        
        configureWeigths(self.canvas)
                
        self.contentFrame = tk.Frame(self.canvas, bg='white', bd=0)
        self.contentId = self.canvas.create_window(0, 0, anchor=tk.NW, window=self.contentFrame)
        
        self.contentFrame.bind('<Configure>', self._configure_interior)
        self.canvas.bind('<Configure>', self._configure_canvas)
        
        self.contentFrame.columnconfigure(0, weight=1)
        #configureWeigths(self.contentFrame)
        self.columnconfigure(0, weight=1)
        
        
    def _getReqSize(self, widget):
        return widget.winfo_reqwidth(), widget.winfo_reqheight()
    
    def _getSize(self, widget):
        return widget.winfo_width(), widget.winfo_height()
    # track changes to the canvas and frame width and sync them,
    # also updating the scrollbar
    def _configure_interior(self, event=None):
        # update the scrollbars to match the size of the inner frame
        fsize = self._getReqSize(self.contentFrame)
        csize = self._getSize(self.canvas)
        #if fsize[0] != self.canvas.winfo_width():
        if fsize != csize:
            # update the canvas's width to fit the inner frame
            self.canvas.config(width=fsize[0], height=fsize[1])
            self.canvas.config(scrollregion="0 0 %s %s" % fsize)

    def _configure_canvas(self, event=None):
        fsize = self._getReqSize(self.contentFrame)
        csize = self._getSize(self.canvas)
        #if self.contentFrame.winfo_reqwidth() != self.canvas.winfo_width():
        if fsize != csize:
            # update the inner frame's width to fill the canvas
            self.canvas.itemconfigure(self.contentId, width=csize[0])
            if csize[1] > fsize[1]:
                self.canvas.itemconfigure(self.contentId, height=csize[1])
                self.canvas.config(scrollregion="0 0 %s %s" % csize)
                
    def adjustContent(self):
        self._configure_interior()
        self.update_idletasks()
        self._configure_canvas()
        
                    
class SectionWidget(SectionFrame):
    """This class will be used to create a section in FormWindow"""
    def __init__(self, form, master, section, height, callback=None, **args):
        self.form = form
        self.section = section
        self.callback = callback
        SectionFrame.__init__(self, master, self.section.label.get(), height=height, **args)
        
    def _createHeader(self, label, bgColor):
        SectionFrame._createHeader(self, label, bgColor)        
        
        if self.section.hasQuestion():
            question = self.section.getQuestion() 
            self.paramName = self.section.getQuestionName()            
            self.var = BoolVar()
            self.var.set(question.get())
            self.var.trace('w', self._onVarChanged)
            
            self.chbLabel = tk.Label(self.headerFrame, text=question.label.get(), fg='white', bg=bgColor)
            self.chbLabel.grid(row=0, column=1, sticky='e', padx=2)
            
            self.chb = tk.Checkbutton(self.headerFrame, variable=self.var.tkVar, 
                                      bg=bgColor, activebackground=gui.cfgButtonActiveBgColor)
            self.chb.grid(row=0, column=2, sticky='e')
        
    def show(self):
        self.contentFrame.grid(row=1, column=0, sticky='news', padx=5, pady=5)
        
    def hide(self):
        self.contentFrame.grid_remove()
        
    def _onVarChanged(self, *args):
        if self.get():
            self.show()
        else:
            self.hide()
            
        if self.callback is not None:
            self.callback(self.paramName) 
            
    def get(self):
        """Return boolean value if is selected"""
        return self.var.get()
    
    def set(self, value):
        self.var.set(value)
    
               
class ParamWidget():
    """For each one in the Protocol parameters, there will be
    one of this in the Form GUI.
    It is mainly composed by:
    A Label: put in the left column
    A Frame(content): in the middle column and container
      of the specific components for this parameter
    A Frame(buttons): a container for available actions buttons
    It will also have a Variable that should be set when creating 
      the specific components"""
    def __init__(self, row, paramName, param, window, parent, value, 
                 callback=None, visualizeCallback=None, column=0, showButtons=True):
        self.window = window
        self._protocol = self.window.protocol
        self.row = row
        self.column = column
        self.paramName = paramName
        self.param = param
        self.parent = parent
        self.visualizeCallback = visualizeCallback
        self.var = None
        
        self._btnCol = 0
        self._labelFont = self.window.font

        self._initialize(showButtons)
        self._createLabel() # self.label should be set after this 
        self._createContent() # self.content and self.var should be set after this
        
        if self.var: # Groups have not self.var
            self.set(value)
            self.callback = callback
            self.var.trace('w', self._onVarChanged)
        
    def _initialize(self, showButtons):
        # Show buttons = False means the widget is inside a Line group
        # then, some of the properties change accordingly
        if showButtons: 
            self._labelSticky = 'ne'
            self._padx, self._pady = 2, 2
            self._entryWidth = 10
            if self.param.isImportant():
                self._labelFont = self.window.fontBold
            self.parent.columnconfigure(0, minsize=250)
            self.parent.columnconfigure(1, minsize=250)
            self.btnFrame = tk.Frame(self.parent, bg='white') # self.btnFrame should be set after this
        else:
            self.btnFrame = None
            self._labelSticky = 'nw'
            self._padx, self._pady = 2, 0
            self._labelFont = self.window.fontItalic
            self._entryWidth = 8
        
    def _createLabel(self):
        bgColor = 'white'
        
        if self.param.isExpert():
            bgColor = 'grey'
        
        self.label = tk.Label(self.parent, text=self.param.label.get(), 
                              bg=bgColor, font=self._labelFont, wraplength=300)
               
    def _createContent(self):
        self.content = tk.Frame(self.parent, bg='white')
        gui.configureWeigths(self.content) 
        self._createContentWidgets(self.param, self.content) # self.var should be set after this
        
    def _addButton(self, text, imgPath, cmd):
        if self.btnFrame:
            btn = IconButton(self.btnFrame, text, imgPath, command=cmd)
            btn.grid(row=0, column=self._btnCol, sticky='e', padx=2, pady=2)
            self.btnFrame.columnconfigure(self._btnCol, weight=1)
            self._btnCol += 1
        
    def _showHelpMessage(self, e=None):
        showInfo("Help", self.param.help.get(), self.parent)
        
    def _showInfo(self, msg):
        showInfo("Info", msg, self.parent)
        
    def _showWizard(self, e=None):
        wizClass = self.window.wizards[self.wizParamName]
        wizard = wizClass()
        wizard.show(self.window)
        
    def _findParamWizard(self):
        """ Search if there are registered wizards for this param
        or any of its subparams (for the case of Line groups)
        """
        if self.paramName in self.window.wizards:
            self.wizParamName = self.paramName
            return True
        
        if isinstance(self.param, Line):
            for name, _ in self.param.iterParams():
                if name in self.window.wizards:
                    self.wizParamName = name
                    return True
        # Search in sub-params
        return False
               
    @staticmethod
    def createBoolWidget(parent, **args):
        """ Return a BoolVar associated with a yes/no selection. 
        **args: extra arguments passed to tk.Radiobutton and tk.Frame constructors.
        """
        var = BoolVar()
        frame = tk.Frame(parent, **args)
        rb1 = tk.Radiobutton(frame, text='Yes', variable=var.tkVar, value=1, **args)
        rb1.grid(row=0, column=0, padx=2, sticky='w')
        rb2 = tk.Radiobutton(frame, text='No', variable=var.tkVar, value=0, **args)
        rb2.grid(row=0, column=1, padx=2, sticky='w') 
        
        return (var, frame)
    
    def _createContentWidgets(self, param, content):
        """Create the specific widgets inside the content frame"""
        # Create widgets for each type of param
        t = type(param)
        entryWidth = 30

        if t is BooleanParam:
            var, frame = ParamWidget.createBoolWidget(content, bg='white')
            frame.grid(row=0, column=0, sticky='w')
            
        elif t is EnumParam:
            var = ComboVar(param)
            if param.display == EnumParam.DISPLAY_COMBO:
                combo = ttk.Combobox(content, textvariable=var.tkVar, state='readonly')
                combo['values'] = param.choices
                combo.grid(row=0, column=0, sticky='w')
            elif param.display == EnumParam.DISPLAY_LIST:
                for i, opt in enumerate(param.choices):
                    rb = tk.Radiobutton(content, text=opt, variable=var.tkVar, value=opt)
                    rb.grid(row=i, column=0, sticky='w')
            else:
                raise Exception("Invalid display value '%s' for EnumParam" % str(param.display))
        
        elif t is MultiPointerParam:
            tp = MultiPointerTreeProvider()
            tree = BoundTree(content, tp)
            var = MultiPointerVar(tp, tree)
            tree.grid(row=0, column=0, sticky='w')
            self._addButton("Select", Icon.ACTION_SEARCH, self._browseObject)
            self._addButton("Remove", Icon.ACTION_DELETE, self._removeObject)
        
        elif t is PointerParam or t is RelationParam:
            var = PointerVar(self._protocol)
            entry = tk.Entry(content, width=entryWidth, textvariable=var.tkVar, state="readonly")
            entry.grid(row=0, column=0, sticky='w')
            btnFunc = self._browseObject
            if t is RelationParam:
                btnFunc = self._browseRelation
            self._addButton("Select", Icon.ACTION_SEARCH, btnFunc)
        
        elif t is ProtocolClassParam:
            var = tk.StringVar()
            entry = tk.Entry(content, width=entryWidth, textvariable=var, state="readonly")
            entry.grid(row=0, column=0, sticky='w')

            protClassName = self.param.protocolClassName.get()
            
            if self.param.allowSubclasses:
                from pyworkflow.em import findSubClasses, emProtocolsDict
                classes = findSubClasses(emProtocolsDict, protClassName).keys()
            else:
                classes = [protClassName]
            
            if len(classes) > 1:
                self._addButton("Select", Icon.ACTION_SEARCH, self._browseProtocolClass)
            else:
                var.set(classes[0])
            
            self._addButton("Edit", Icon.ACTION_EDIT, self._openProtocolForm)
            #btn = Button(content, "Edit", command=self._openProtocolForm)
            #btn.grid(row=1, column=0)
        elif t is Line:
            var = None
        else:
            #v = self.setVarValue(paramName)
            var = tk.StringVar()
            if issubclass(t, FloatParam) or issubclass(t, IntParam):
                entryWidth = self._entryWidth # Reduce the entry width for numbers entries
            entry = tk.Entry(content, width=entryWidth, textvariable=var)
            entry.grid(row=0, column=0, sticky='w')
            
            if issubclass(t, PathParam):
                self._addButton('Browse', Icon.ACTION_BROWSE, self._browsePath)

        if self.visualizeCallback is not None:
            self._addButton(Message.LABEL_BUTTON_VIS, Icon.ACTION_VISUALIZE, self._visualizeVar)    
        
        if self._findParamWizard():
            self._addButton(Message.LABEL_BUTTON_WIZ, Icon.ACTION_WIZ, self._showWizard)
        
        if param.help.hasValue():
            self._addButton(Message.LABEL_BUTTON_HELP, Icon.ACTION_HELP, self._showHelpMessage)
        
        self.var = var
        
    def _visualizeVar(self, e=None):
        """ Visualize specific variable. """
        self.visualizeCallback(self.paramName)
        
    def _browseObject(self, e=None):
        """Select an object from DB
        This function is suppose to be used only for PointerParam"""
        value = self.get()
        selected = []
        if isinstance(value, list):
            selected = value
        elif selected is not None:
            selected = [value]
        tp = SubclassesTreeProvider(self._protocol, self.param, selected=selected)
        dlg = ListDialog(self.parent, "Select object", tp, 
                         "Double click an item to preview the object")
        if dlg.value is not None:
            self.set(dlg.value)
        
    def _removeObject(self, e=None):
        """ Remove an object from a MultiPointer param. """
        self.var.remove()
                        
    def _browseRelation(self, e=None):
        """Select a relation from DB
        This function is suppose to be used only for RelationParam. """
        tp = RelationsTreeProvider(self._protocol, self.param, selected=self.get())
        dlg = ListDialog(self.parent, "Select object", tp)
        if dlg.value is not None:
            self.set(dlg.value)
            
    def _browseProtocolClass(self, e=None):
        tp = ProtocolClassTreeProvider(self.param.protocolClassName.get())
        dlg = ListDialog(self.parent, "Select protocol", tp)
        if dlg.value is not None:
            self.set(dlg.value)
            self._openProtocolForm()
            
    def _browsePath(self, e=None):
        def onSelect(obj):
            self.set(obj.getPath())
        v = self.get().strip()
        path = None
        if v:
            v = os.path.dirname(v)
            if os.path.exists(v):
                path = v        
        if not path:
            path = getHomePath()
        browser = FileBrowserWindow("Browsing", self.window, path=path, onSelect=onSelect)
        browser.show()
            
    def _openProtocolForm(self, e=None):
        className = self.get().strip()
        
        if len(className):
            instanceName = self.paramName + "Instance"
            protocol = self._protocol
            #TODO check if is present and is selected a different
            # class, so we need to delete that and create a new instance
            if not hasattr(protocol, instanceName):
                from pyworkflow.em import findClass
                cls = findClass(className)
                protocol._insertChild(instanceName, cls())
            
            prot = getattr(protocol, instanceName)
                
            prot.allowHeader.set(False)
            f = FormWindow("Sub-Protocol: " + instanceName, prot, self._protocolFormCallback, self.window, childMode=True)
            f.show()
        else:
            self._showInfo("Select the protocol class first")
        
    def _protocolFormCallback(self, e=None):
        pass
    
    def _onVarChanged(self, *args):
        if self.callback is not None:
            self.callback(self.paramName)        
        
    def show(self):
        """Grid the label and content in the specified row"""
        c = self.column
        self.label.grid(row=self.row, column=c, sticky=self._labelSticky, padx=self._padx, pady=self._pady)
        self.content.grid(row=self.row, column=c+1, sticky='news', 
                          padx=self._padx, pady=self._pady)
        if self.btnFrame:
            self.btnFrame.grid(row=self.row, column=c+2, padx=self._padx, sticky='new')
        
    def hide(self):
        self.label.grid_remove()
        self.content.grid_remove()
        if self.btnFrame:
            self.btnFrame.grid_remove()
            
    def display(self, condition):
        """ show or hide depending on the condition. """
        if condition:
            self.show()
        else:
            self.hide()
        
    def set(self, value):
        if value is not None:
            self.var.set(value)
        
    def get(self):
        return self.var.get()


class LineWidget(ParamWidget):
    def __init__(self, row, paramName, param, window, parent, value, 
                 callback=None, visualizeCallback=None, column=0, showButtons=True):
        ParamWidget.__init__(self, row, paramName, param, window, parent, None)
        self.show()
        
    def show(self):
        self.label.grid(row=self.row, column=0, sticky=self._labelSticky)
        self.content.grid(row=self.row, column=1, sticky='nw', columnspan=6, padx=5)
        if self.btnFrame:
            self.btnFrame.grid(row=self.row, column=2, padx=2, sticky='new')
       

class GroupWidget(ParamWidget):
    def __init__(self, row, paramName, param, window, parent):
        ParamWidget.__init__(self, row, paramName, param, window, parent, None)
        
    def _initialize(self, showButtons):
        pass
        
    def _createLabel(self):
        pass
               
    def _createContent(self):
        self.content = tk.LabelFrame(self.parent, text=self.paramName, bg='white')
        gui.configureWeigths(self.content) 
        
    def show(self):
        self.content.grid(row=self.row, column=0, sticky='news', columnspan=6, padx=5, pady=5)
        
    def hide(self):
        self.content.grid_remove()  
            
            
class Binding():
    def __init__(self, paramName, var, protocol, *callbacks):
        self.paramName = paramName
        self.var = var
        self.var.set(protocol.getAttributeValue(paramName, ''))
        self.var.trace('w', self._onVarChanged)
        self.callbacks = callbacks
        
    def _onVarChanged(self, *args):
        for cb in self.callbacks:
            cb(self.paramName)
            
    
class FormWindow(Window):
    """ This class will create the Protocol params GUI to fill in the parameters.
    The creation of input parameters will be based on the Protocol Form definition.
    This class will serve as a connection between the GUI variables (tk vars) and 
    the Protocol variables.
    
    Layout:
        There are 4 main blocks that goes each one in a different row1.
        1. Header: will contains the logo, title and some link buttons.
        2. Common: common execution parameters of each run.
        3. Params: the expert level and tabs with the Protocol parameters.
        4. Buttons: buttons at bottom for close, save and execute.
    """
    def __init__(self, title, protocol, callback, master=None, hostList=['localhost'], **args):
        """ Constructor of the Form window. 
        Params:
         title: title string of the windows.
         protocol: protocol from which the form will be generated.
         callback: callback function to call when Save or Excecute are press.
        """
        Window.__init__(self, title, master, icon='scipion_bn.xbm', 
                        weight=False, minsize=(600, 450), **args)

        # Some initialization
        self.callback = callback
        self.widgetDict = {} # Store tkVars associated with params
        self.visualizeDict = args.get('visualizeDict', {})
        self.bindings = []
        self.hostList = hostList
        self.protocol = protocol
        self.visualizeMode = args.get('visualizeMode', False)  # This control when to close or not after execute
        self.headerBgColor = Color.RED_COLOR
        if self.visualizeMode:
            self.headerBgColor = Color.DARK_GREY_COLOR
        self.childMode = args.get('childMode', False) # Allow to open child protocols form (for workflows)
        
        from pyworkflow.em import findWizards
        self.wizards = findWizards(protocol, DESKTOP_TKINTER)
        
        self._createGUI()
        
    def _createGUI(self):
        mainFrame = tk.Frame(self.root)
        configureWeigths(mainFrame, row=2)
        self.root.rowconfigure(0, weight=1)
        
        headerFrame = self._createHeader(mainFrame)
        headerFrame.grid(row=0, column=0, sticky='new')
        
        if self.protocol.allowHeader:
            commonFrame = self._createCommon(mainFrame)
            commonFrame.grid(row=1, column=0, sticky='nw')
            
        paramsFrame = self._createParams(mainFrame)
        paramsFrame.grid(row=2, column=0, sticky='news')
        
        buttonsFrame = self._createButtons(mainFrame)
        buttonsFrame.grid(row=3, column=0, sticky='se')
        
        mainFrame.grid(row=0, column=0, sticky='ns')
        
        
    def _createHeader(self, parent):
        """ Fill the header frame with the logo, title and cite-help buttons. """
        headerFrame = tk.Frame(parent)
        #headerFrame.grid(row=0, column=0, sticky='new')
        headerFrame.columnconfigure(0, weight=1)
        package = self.protocol._package
        t = '  Protocol: %s' % (self.protocol.getClassLabel())
        logoPath = getattr(package, '_logo', 'scipion_bn.png')
        
        if logoPath:
            headerLabel = tk.Label(headerFrame, text=t, font=self.fontBig, 
                                   image=self.getImage(logoPath, maxheight=40), compound=tk.LEFT)
        else:
            headerLabel = tk.Label(headerFrame, text=t, font=self.fontBig)
        headerLabel.grid(row=0, column=0, padx=5, pady=(5,0), sticky='nw')#, columnspan=5)
        
        def _addButton(text, icon, command, col):
            btn = tk.Label(headerFrame, text=text, image=self.getImage(icon), 
                       compound=tk.LEFT, cursor='hand2')
            btn.bind('<Button-1>', command)
            btn.grid(row=0, column=col, padx=5, sticky='e')
        
        _addButton(Message.LABEL_CITE, Icon.ACTION_REFERENCES, self._showReferences, 1)
        _addButton(Message.LABEL_HELP ,Icon.ACTION_HELP, self._showHelp, 2)
        
        return headerFrame
        
    def _showReferences(self, e=None):
        """ Show the list of references of the protocol. """
        self.showInfo('\n'.join(self.protocol.citations()), "References")
        
    def _showHelp(self, e=None):
        """ Show the list of references of the protocol. """
        self.showInfo(self.protocol.getDoc(), "Help")
        
        
    def _createCommon(self, parent):
        """ Create the second section with some common parameters. """
        commonFrame = tk.Frame(parent)
        
        ############# Create the run part ###############
        runSection = SectionFrame(commonFrame, label=Message.TITLE_RUN, height=100,
                                  headerBgColor=self.headerBgColor)
        runFrame = tk.Frame(runSection.contentFrame, bg='white')
        runFrame.grid(row=0, column=0, sticky='nw', padx=(0, 15))
        #runFrame.columnconfigure(1, weight=1)
        
        r = 0 # Run name
        self._createHeaderLabel(runFrame, Message.LABEL_RUNNAME, bold=True, sticky='ne')
        self.runNameVar = tk.StringVar()
        entry = tk.Entry(runFrame, font=self.font, width=25, textvariable=self.runNameVar)
        entry.grid(row=r, column=1, padx=(0, 5), pady=5, sticky='new')#, columnspan=5)
        btn = IconButton(runFrame, Message.TITLE_COMMENT, Icon.ACTION_EDIT, command=self._editObjParams)
        btn.grid(row=r, column=2, padx=(10,0), pady=5, sticky='nw')
        
        c = 3 # Comment
        self._createHeaderLabel(runFrame, Message.TITLE_COMMENT, sticky='ne', column=c)
        self.commentVar = tk.StringVar()
        entry = tk.Entry(runFrame, font=self.font, width=25, textvariable=self.commentVar)
        entry.grid(row=r, column=c+1, padx=(0, 5), pady=5, sticky='new')#, columnspan=5)
        btn = IconButton(runFrame, Message.TITLE_COMMENT, Icon.ACTION_EDIT, command=self._editObjParams)
        btn.grid(row=r, column=c+2, padx=(10,0), pady=5, sticky='nw')
        
        self.updateLabelAndComment()
                
        r = 1 # Execution
        self._createHeaderLabel(runFrame, Message.LABEL_EXECUTION, bold=True, sticky='ne', row=r, pady=0)
        modeFrame = tk.Frame(runFrame, bg='white')
        self._createHeaderLabel(modeFrame, "Mode", sticky='ne', row=0, pady=0, column=0)
        runMode = self._createBoundCombo(modeFrame, Message.VAR_RUN_MODE, MODE_CHOICES, self._onRunModeChanged)   
        runMode.grid(row=0, column=1, sticky='new', padx=(0, 5), pady=5)
        btnHelp = IconButton(modeFrame, Message.TITLE_COMMENT, Icon.ACTION_HELP, 
                             command=self._createHelpCommand(Message.HELP_RUNMODE))
        btnHelp.grid(row=0, column=2, padx=(5, 0), pady=2, sticky='ne')
        modeFrame.columnconfigure(0, minsize=60)
        modeFrame.columnconfigure(1, weight=1)
        modeFrame.grid(row=r, column=1, sticky='new', columnspan=2)
        
        # Host
        self._createHeaderLabel(runFrame, Message.LABEL_HOST, row=r, column=c, pady=0, padx=(15,5), sticky='ne')
        self._createBoundCombo(runFrame, Message.VAR_EXEC_HOST, self.hostList).grid(row=r, column=c+1, pady=5, sticky='nw')
        
        r = 2
        # Parallel
        if self.protocol.allowThreads or self.protocol.allowMpi:
            self._createHeaderLabel(runFrame, Message.LABEL_PARALLEL, bold=True, sticky='ne', row=r, pady=0)
            procFrame = tk.Frame(runFrame, bg='white')
            r2 = 0
            c2 = 0
            sticky = 'ne'
            # THREADS
            if self.protocol.allowThreads:
                self._createHeaderLabel(procFrame, Message.LABEL_THREADS, sticky=sticky, row=r2, column=c2, pady=0)
                entry = self._createBoundEntry(procFrame, Message.VAR_THREADS)
                entry.grid(row=r2, column=c2+1, padx=(0, 5), sticky='nw')
                # Modify values to be used in MPI entry
                c2 += 2
                sticky = 'nw'
            # MPI
            if self.protocol.allowMpi:
                self._createHeaderLabel(procFrame, Message.LABEL_MPI, sticky=sticky, row=r2, column=c2, pady=0)
                entry = self._createBoundEntry(procFrame, Message.VAR_MPI)
                entry.grid(row=r2, column=c2+1, padx=(0, 5), sticky='nw')
                
            btnHelp = IconButton(procFrame, Message.TITLE_COMMENT, Icon.ACTION_HELP, 
                                 command=self._createHelpCommand(Message.HELP_MPI_THREADS))
            btnHelp.grid(row=0, column=4, padx=(5, 0), pady=2, sticky='ne')
            procFrame.columnconfigure(0, minsize=60)
            procFrame.grid(row=r, column=1, sticky='new', columnspan=2)
        
        # Queue
        self._createHeaderLabel(runFrame, Message.LABEL_QUEUE, row=r, sticky='ne', column=c, padx=(15,5), pady=0)
        var, frame = ParamWidget.createBoolWidget(runFrame, bg='white')
        self._addVarBinding(Message.VAR_QUEUE, var)
        frame.grid(row=2, column=c+1, pady=5, sticky='nw')
        btnHelp = IconButton(runFrame, Message.TITLE_COMMENT, Icon.ACTION_HELP, 
                             command=self._createHelpCommand(Message.HELP_RUNMODE))
        btnHelp.grid(row=2, column=c+2, padx=(5, 0), pady=2, sticky='ne')
        
        # Run Name not editable
        #entry.configure(state='readonly')
        # Run mode
        self.protocol.getDefinitionParam('')
        #self._createHeaderLabel(runFrame, Message.LABEL_RUNMODE).grid(row=1, column=0, sticky='ne', padx=5, pady=5)
        #runSection.addContent()
        runSection.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        
        return commonFrame 
 
    def _createHelpCommand(self, msg):
        """ Show the help of some value of the header. """
        return lambda: showInfo("Help", msg, self.root)
    
    def _editObjParams(self, e=None):
        """ Show a Text area to edit the protocol label and comment. """
        self.setProtocolLabel()        
        d = EditObjectDialog(self.root, Message.TITLE_EDIT_OBJECT, self.protocol, self.protocol.mapper)
        
        if d.resultYes():
            self.updateLabelAndComment()
        
    def _createParams(self, parent):
        paramsFrame = tk.Frame(parent)
        configureWeigths(paramsFrame, row=1, column=0)
        # Expert level
        expFrame = tk.Frame(paramsFrame)
        expLabel = tk.Label(expFrame, text=Message.LABEL_EXPERT, font=self.fontBold)
        expLabel.grid(row=0, column=0, sticky='nw', padx=5)
        expCombo = self._createBoundCombo(expFrame, Message.VAR_EXPERT, LEVEL_CHOICES, self._onExpertLevelChanged) 
        expCombo.grid(row=0, column=1, sticky='nw', pady=5)
        expFrame.grid(row=0, column=0, sticky='nw')
        
        contentFrame = self._createSections(paramsFrame)
        contentFrame.grid(row=1, column=0, sticky='news')
        
        return paramsFrame
                
    def _createSections(self, parent):
        """Create section widgets"""
        r = 0
        sectionsFrame = tk.Frame(parent) 
        configureWeigths(sectionsFrame)
        tab = ttk.Notebook(sectionsFrame) 
        tab.grid(row=0, column=0, sticky='news',
                 padx=5, pady=5)
        self._sections = []
        
        for section in self.protocol.iterDefinitionSections():
            label = section.getLabel()
            if label != 'General' and label != 'Parallelization':
                frame = SectionWidget(self, tab, section, height=150,
                                      callback=self._checkChanges,
                                      headerBgColor=self.headerBgColor)
                
                tab.add(frame, text=section.getLabel())
                frame.columnconfigure(0, minsize=400)
                self._fillSection(section, frame)
                self._sections.append(frame)
                r += 1
        self._checkAllChanges()
        
        return sectionsFrame    
        
    def _createButtons(self, parent):
        """ Create the bottom buttons: Close, Save and Execute. """
        btnFrame = tk.Frame(parent)
        
        btnClose = self.createCloseButton(btnFrame)
        btnClose.grid(row=0, column=0, padx=5, pady=5, sticky='se')
        # Save button is not added in VISUALIZE or CHILD modes
        if not self.visualizeMode and not self.childMode:
            btnSave = Button(btnFrame, Message.LABEL_BUTTON_RETURN, Icon.ACTION_SAVE, 
                                command=self.save)
            btnSave.grid(row=0, column=1, padx=5, pady=5, sticky='se')
            btnExecute = HotButton(btnFrame, Message.LABEL_BUTTON_EXEC, 
                                   Icon.ACTION_EXECUTE, command=self.execute)
            btnExecute.grid(row=0, column=2, padx=(5, 28), pady=5, sticky='se')
            
        return btnFrame
        
    def _addVarBinding(self, paramName, var, func=None, *callbacks):
        if func is None:
            func = self.setParamFromVar
        binding = Binding(paramName, var, self.protocol, 
                          func, *callbacks)
        self.widgetDict[paramName] = var
        self.bindings.append(binding)
        
    def _createBoundEntry(self, parent, paramName, width=5, func=None, value=None):
        var = tk.StringVar()
        setattr(self, paramName + 'Var', var)
        self._addVarBinding(paramName, var, func)
        if value is not None:
            var.set(value)
        return tk.Entry(parent, font=self.font, width=width, textvariable=var)
    
    def _createEnumBinding(self, paramName, choices, *callbacks):
        param = EnumParam(choices=choices)
        var = ComboVar(param)
        self._addVarBinding(paramName, var, None, *callbacks)
        return param, var
        
    def _createBoundCombo(self, parent, paramName, choices, *callbacks):
        param, var = self._createEnumBinding(paramName, choices, *callbacks)
        combo = ttk.Combobox(parent, textvariable=var.tkVar, state='readonly', width=10)
        combo['values'] = param.choices
        
        return combo
    
    def _createBoundOptions(self, parent, paramName, choices, *callbacks, **args):
        param, var = self._createEnumBinding(paramName, choices, *callbacks)
        frame = tk.Frame(parent, **args)
        for i, opt in enumerate(param.choices):
            rb = tk.Radiobutton(frame, text=opt, variable=var, value=opt)
            rb.grid(row=0, column=i, sticky='nw', padx=(0, 5))  
        
        return frame
        
    def _createHeaderLabel(self, parent, text, bold=False, **gridArgs):
        font = self.font
        if bold:
            font = self.fontBold
        label = tk.Label(parent, text=text, font=font, bg='white')
        if gridArgs:
            gridDefaults = {'row': 0, 'column': 0, 'padx': 5, 'pady': 5}
            gridDefaults.update(gridArgs)
            label.grid(**gridDefaults)
        return label
    
    def resize(self, frame):
        self.root.update_idletasks()
        MaxHeight = 1200
        MaxWidth = 1600
        rh = frame.winfo_reqheight()
        rw = frame.winfo_reqwidth()
        height = min(rh + 100, MaxHeight)
        width = min(rw, MaxWidth)
        x = self.root.winfo_x()
        y = self.root.winfo_y()
        self.root.geometry("%dx%d%+d%+d" % (width, height, x, y))

        return (width, height)
    
    def adjustSize(self):
        self.resize(self.root)        
        
    def save(self, e=None):
        self._close(onlySave=True)
        
    def execute(self, e=None):
        if (self.protocol.getRunMode() == MODE_RESTART and 
            not askYesNo(Message.TITLE_RESTART_FORM, Message.LABEL_DELETE_FORM, self.root)):
            return
            
        errors = self.protocol.validate()
        
        if len(errors):
            self.showError(errors)
        else:
            self._close()
        
    def _close(self, onlySave=False):
        try:
            # Set the protocol label
            self.setProtocolLabel()
            
            message = self.callback(self.protocol, onlySave)
            if not self.visualizeMode:
                if len(message):
                    self.showInfo(message, "Protocol action")
                if not onlySave:
                    self.close()
        except Exception, ex:
            action = "EXECUTE"
            if onlySave:
                action = "SAVE"
            self.showError("Error during %s: %s" % (action, ex))
    
    
    def getWidgetValue(self, protVar, param):
        widgetValue = ""                
        if isinstance(param, MultiPointerParam):
            for pointer in protVar:
                widgetValue += pointer.get(param.default.get()) + ";"
        else:
            widgetValue = protVar.get(param.default.get())  
        return widgetValue
          
    def _visualize(self, paramName):
        protVar = getattr(self.protocol, paramName)
        if protVar.hasValue():
            from pyworkflow.em import findViewers
            obj = protVar.get() # Get the reference to the object
            viewers = findViewers(obj.getClassName(), DESKTOP_TKINTER)
            if len(viewers):
                v = viewers[0] # Use the first viewer registered
                proj = self.protocol.getProject()
                if proj is None:
                    raise Exception("none project")
                v(project=self.protocol.getProject()).visualize(obj) # Instanciate the viewer and visualize object
            else:
                self.showInfo("There is not viewer registered for this object")
        else:
            self.showInfo("Select the object before visualize")
         
    def _fillSection(self, sectionParam, sectionWidget):
        parent = sectionWidget.contentFrame
        r = 0
        for paramName, param in sectionParam.iterParams():
            if isinstance(param, Group):
                widget = GroupWidget(r, paramName, param, self, parent)
                self._fillGroup(param, widget)
            elif isinstance (param, Line):
                widget = LineWidget(r, paramName, param, self, parent, None)
                self._fillLine(param, widget)
            else:
                protVar = getattr(self.protocol, paramName, None)
                
                if protVar is None:
                    raise Exception("_fillSection: param '%s' not found in protocol" % paramName)
                
                if sectionParam.getQuestionName() == paramName:
                    widget = sectionWidget
                    if not protVar:
                        widget.hide() # Show only if question var is True
                else:
                    if isinstance(param, PointerParam):
                        visualizeCallback = self._visualize # Add visualize icon for pointer params
                    else:
                        visualizeCallback = self.visualizeDict.get(paramName, None)
                    
                    widget = ParamWidget(r, paramName, param, self, parent, 
                                                             value=self.getWidgetValue(protVar, param),
                                                             callback=self._checkChanges,
                                                             visualizeCallback=visualizeCallback)
                        
                    widget.show() # Show always, conditions will be checked later
            r += 1         
            self.widgetDict[paramName] = widget
        # Ensure width and height needed
        w, h = parent.winfo_reqwidth(), parent.winfo_reqheight()
        sectionWidget.columnconfigure(0, minsize=w)
        sectionWidget.rowconfigure(0, minsize=h)

    def _fillGroup(self, groupParam, groupWidget):
        parent = groupWidget.content
        r = 0
        for paramName, param in groupParam.iterParams():
            protVar = getattr(self.protocol, paramName, None)
            
            if protVar is None:
                raise Exception("_fillSection: param '%s' not found in protocol" % paramName)
            
            if isinstance(param, PointerParam):
                visualizeCallback = self._visualize # Add visualize icon for pointer params
            else:
                visualizeCallback = self.visualizeDict.get(paramName, None)
            
            widget = ParamWidget(r, paramName, param, self, parent, 
                                                     value=self.getWidgetValue(protVar, param),
                                                     callback=self._checkChanges,
                                                     visualizeCallback=visualizeCallback)
            widget.show() # Show always, conditions will be checked later
            r += 1         
            self.widgetDict[paramName] = widget
 
    def _fillLine(self, groupParam, groupWidget):
        parent = groupWidget.content
        r = 0
        for paramName, param in groupParam.iterParams():
            protVar = getattr(self.protocol, paramName, None)
            
            if protVar is None:
                raise Exception("_fillSection: param '%s' not found in protocol" % paramName)
            
            if isinstance(param, PointerParam):
                visualizeCallback = self._visualize # Add visualize icon for pointer params
            else:
                visualizeCallback = self.visualizeDict.get(paramName, None)
            
            widget = ParamWidget(0, paramName, param, self, parent, 
                                 value=self.getWidgetValue(protVar, param),
                                 callback=self._checkChanges, visualizeCallback=visualizeCallback,
                                 column=r, showButtons=False)
            widget.show() # Show always, conditions will be checked later
            r += 2         
            self.widgetDict[paramName] = widget           

        
    def _checkCondition(self, paramName):
        """Check if the condition of a param is statisfied 
        hide or show it depending on the result"""
        widget = self.widgetDict.get(paramName, None)
        
        if isinstance(widget, ParamWidget): # Special vars like MPI, threads or runName are not real widgets
            if isinstance(widget, LineWidget) or isinstance(widget, GroupWidget):
                param = widget.param
            else:
                param = self.protocol.getDefinitionParam(paramName)
            cond = self.protocol.evalParamCondition(paramName) and self.protocol.evalParamExpertLevel(param)
            widget.display(cond)
            
    def _checkChanges(self, paramName):
        """Check the conditions of all params affected
        by this param"""
        self.setParamFromVar(paramName)
        param = self.protocol.getDefinitionParam(paramName)
        
        for d in param._dependants:
            self._checkCondition(d)
            
    def _checkAllChanges(self):
        for paramName in self.widgetDict:
            self._checkCondition(paramName)
            
    def _onExpertLevelChanged(self, *args):
        self._checkAllChanges()
        self.root.update_idletasks()
        for s in self._sections:
            s.adjustContent()
        
    def _onRunModeChanged(self, paramName):
        self.setParamFromVar(paramName)
        
    def getVarValue(self, varName):
        """This method should retrieve a value from """
        pass
        
    def setVar(self, paramName, value):
        var = self.widgetDict[paramName]
        var.set(value)
        
    def setVarFromParam(self, paramName):
        var = self.widgetDict[paramName]
        param = getattr(self.protocol, paramName, None)
        if param is not None:
            var.set(param.get(''))
           
    def setParamFromVar(self, paramName):
        param = getattr(self.protocol, paramName, None)
        if param is not None:
            var = self.widgetDict[paramName]
            try:
                param.set(var.get())
            except ValueError:
                if len(var.get()):
                    print "Error setting param for: ", paramName, "value: '%s'" % var.get()
                param.set(None)
                
    def updateLabelAndComment(self):
        """ Read the label and comment first line to update
        the entry boxes in the form.
        """
        self.runNameVar.set(self.protocol.getObjLabel())
        # Get only the first comment line
        comment = self.protocol.getObjComment()
        if comment:
            lines = comment.split('\n')
            if lines:
                comment = lines[0]
        self.commentVar.set(comment)
        
    def setProtocolLabel(self):
        self.protocol.setObjLabel(self.runNameVar.get())
             
    def updateProtocolParams(self):
        for paramName, _ in self.protocol.iterDefinitionAttributes():
            self.setParamFromVar(paramName)


def editObject(self, title, root, obj, mapper):
    """ Show a Text area to edit the protocol label and comment. """    
    return EditObjectDialog(root, title, obj, mapper)
    

if __name__ == '__main__':
    # Just for testing
    from pyworkflow.em import ProtImportMicrographs
    p = ProtImportMicrographs()
    p.sphericalAberration.set(2.3)
    p.setSamplingRate('5.4')
    w = FormWindow(p)
    w.show()
    
   


