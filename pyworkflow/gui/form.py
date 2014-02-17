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

from pyworkflow.mapper.mapper import Mapper
from pyworkflow.mapper import mapper

import Tkinter as tk
import ttk
import tkFont

import gui
from gui import configureWeigths, Window
from text import TaggedText
from widgets import Button
from pyworkflow.protocol.params import *
from pyworkflow.protocol import Protocol
from dialog import showInfo, TextDialog, ListDialog
from tree import TreeProvider
from pyworkflow.utils.properties import Message, Icon, Color
from pyworkflow.viewer import DESKTOP_TKINTER
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
    def __init__(self):
        self.tkVar = tk.StringVar()
        self.value = None
        self.trace = self.tkVar.trace
        
    def set(self, value):
        self.value = value
        v = ''

        label = value.getObjLabel()

        if len(label) > 0:
            v = label
        else:
            v = '%s.%s' % (value.getName(), value.strId())
            
        self.tkVar.set(v)   
            
    def get(self):
        return self.value
    
   
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
        self.tkVar.set(self.enum.choices[value])
                    
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
        objs = []
        for objClass in self.className.split(","):
            for obj in self.mapper.selectByClass(objClass.strip(), objectFilter=self.objFilter):
                objs.append(obj)
        return objs        
#        return self.mapper.selectByClass(self.className, objectFilter=self.objFilter)
            
    def objFilter(self, obj):
        result = True
        # Do not allow to select objects that are childs of the protocol
        if self.protocol.getObjId() == obj.getObjParentId():
            result = False
        # Check that the condition is met
        elif self.condition:
            result = obj.evalCondition(self.condition)
        return result
        
    def getColumns(self):
        return [('Object', 400), ('Info', 250)]
    
    def getObjectInfo(self, obj):
        objName = obj.getNameId()
        selected = self.selected is not None and self.selected.getObjId() == obj.getObjId()
        return {'key': objName, 'values': (str(obj),), 'selected': selected}

    def getObjectActions(self, obj):
        if isinstance(obj, Pointer):
            obj = obj.getName()
        actions = []    
        from pyworkflow.em import findViewers
        viewers = findViewers(obj.getClassName(), DESKTOP_TKINTER)
        for v in viewers:
            actions.append(('Open with %s' % v.__name__, lambda : v().visualize(obj)))
            
        return actions
    
    
#TODO: check if need to inherit from SubclassesTreeProvider
class RelationsTreeProvider(SubclassesTreeProvider):
    """Will implement the methods to provide the object info
    of subclasses objects(of className) found by mapper"""
    def __init__(self, protocol, relationParam):
        self.mapper = protocol.mapper
        parentObject = protocol.getAttributeValue(relationParam.relationParent.get())
        if parentObject is not None:
            queryFunc =  protocol.mapper.getRelationChilds      
            if relationParam.relationReverse:
                queryFunc =  protocol.mapper.getRelationParents
            self.getObjects = lambda: queryFunc(relationParam.relationName.get(), parentObject)
        else:
            self.getObjects = lambda: []   
            
   
#---------------------- Other widgets ----------------------------------------
         
class SectionFrame(tk.Frame):
    """This class will be used to create a frame for the Section
    That will have a header with red color and a content frame
    with white background
    """
    def __init__(self, master, label, callback=None, **args):
        headerBgColor = args.get('headerBgColor', gui.cfgButtonBgColor)
        if 'headerBgColor' in args:
            del args['headerBgColor']
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
        self.contentFrame = tk.Frame(self, bg='white', bd=0)
        self.contentFrame.grid(row=1, column=0, sticky='news', padx=5, pady=5)
        configureWeigths(self.contentFrame)
        self.columnconfigure(0, weight=1)
        
                    
class SectionWidget(SectionFrame):
    """This class will be used to create a section in FormWindow"""
    def __init__(self, form, master, section, callback=None, **args):
        self.form = form
        self.section = section
        self.callback = callback
        SectionFrame.__init__(self, master, self.section.label.get(), **args)
        
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
    def __init__(self, row, paramName, param, window, parent, value, callback=None, visualizeCallback=None):
        self.window = window
        self.row = row
        self.paramName = paramName
        self.param = param
        self.parent = parent
        self.visualizeCallback = visualizeCallback
        
        self.parent.columnconfigure(0, minsize=250)
        self.parent.columnconfigure(1, minsize=250)
        self._createLabel() # self.label should be set after this 
        self._btnCol = 0
        self._createButtonsFrame() # self.btnFrame should be set after this
        self._createContent() # self.content and self.var should be set after this
        
        self.set(value)
        self.callback = callback
        self.var.trace('w', self._onVarChanged)
        
        
    def _createLabel(self):
        f = self.window.font

        if self.param.isImportant:
            f = self.window.fontBold
            
        bgColor = 'white'
        
        if self.param.isExpert():
            bgColor = 'grey'
        
        self.label = tk.Label(self.parent, text=self.param.label.get(), 
                              bg=bgColor, font=f, wraplength=300)
               
    def _createContent(self):
        self.content = tk.Frame(self.parent, bg='white')
        gui.configureWeigths(self.content) 
        self._createContentWidgets(self.param, self.content) # self.var should be set after this
        
    def _addButton(self, text, imgPath, cmd):
        btn = Button(self.btnFrame, text, imgPath, bg='white', command=cmd)
        btn.grid(row=0, column=self._btnCol, sticky='e', padx=2)
        self.btnFrame.columnconfigure(self._btnCol, weight=1)
        self._btnCol += 1
        
    def _createButtonsFrame(self):
        self.btnFrame = tk.Frame(self.parent, bg='white')
        
    def _showHelpMessage(self, e=None):
        showInfo("Help", self.param.help.get(), self.parent)
        
    def _showInfo(self, msg):
        showInfo("Info", msg, self.parent)
        
    def _showWizard(self, e=None):
        wizClass = self.window.wizards[self.paramName]
        wizClass().show(self.window)
               
    @staticmethod
    def createBoolWidget(parent, **args):
        """ Return a BoolVar associated with a yes/no selection. 
        **args: extra arguments passed to tk.Radiobutton and tk.Frame constructors.
        """
        var = BoolVar()
        frame = tk.Frame(parent, **args)
        frame.grid(row=0, column=0, sticky='w')
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
        #TODO: Move this to a Renderer class to be more flexible
        if t is BooleanParam:
            var, frame = ParamWidget.createBoolWidget(content, bg='white')
#            var = BoolVar()
#            frame = tk.Frame(content, bg='white')
#            frame.grid(row=0, column=0, sticky='w')
#            rb1 = tk.Radiobutton(frame, text='Yes', bg='white', variable=var.tkVar, value=1)
#            rb1.grid(row=0, column=0, padx=2, sticky='w')
#            rb2 = tk.Radiobutton(frame, text='No', bg='white', variable=var.tkVar, value=0)
#            rb2.grid(row=0, column=1, padx=2, sticky='w')
            
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
#            listBox = Listbox(content)
#            print "poraki"
            pass
        
        elif t is PointerParam or t is RelationParam:
            var = PointerVar()
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
        else:
            #v = self.setVarValue(paramName)
            var = tk.StringVar()
            if t is FloatParam or t is IntParam:
                entryWidth = 10 # Reduce the entry width for numbers entries
            entry = tk.Entry(content, width=entryWidth, textvariable=var)
            entry.grid(row=0, column=0, sticky='w')

        if self.visualizeCallback is not None:
            self._addButton(Message.LABEL_BUTTON_VIS, Icon.ACTION_VISUALIZE, self._visualizeVar)    
        if self.paramName in self.window.wizards:
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
        tp = SubclassesTreeProvider(self.window.protocol, self.param, selected=self.get())
        dlg = ListDialog(self.parent, "Select object", tp, "Double click an item to preview the object")
        if dlg.value is not None:
            self.set(dlg.value)
            
    def _browseRelation(self, e=None):
        """Select a relation from DB
        This function is suppose to be used only for RelationParam"""
        tp = RelationsTreeProvider(self.window.protocol, self.param)
        dlg = ListDialog(self.parent, "Select object", tp)
        if dlg.value is not None:
            self.set(dlg.value)
            
    def _browseProtocolClass(self, e=None):
        tp = ProtocolClassTreeProvider(self.param.protocolClassName.get())
        dlg = ListDialog(self.parent, "Select protocol", tp)
        if dlg.value is not None:
            self.set(dlg.value)
            self._openProtocolForm()
            
    def _openProtocolForm(self, e=None):
        className = self.get().strip()
        
        if len(className):
            instanceName = self.paramName + "Instance"
            protocol = self.window.protocol
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
        self.label.grid(row=self.row, column=0, sticky='ne', padx=2, pady=2)
        self.content.grid(row=self.row, column=1, padx=2, pady=2, sticky='news')
        self.btnFrame.grid(row=self.row, column=2, padx=2, sticky='new')
        
    def hide(self):
        self.label.grid_remove()
        self.content.grid_remove()
        self.btnFrame.grid_remove()
        
    def set(self, value):
        if value is not None:
            self.var.set(value)
        
    def get(self):
        return self.var.get()
        

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
    """This class will create the Protocol params GUI to enter parameters.
    The creaation of compoments will be based on the Protocol Form definition.
    This class will serve as a connection between the GUI variables (tk vars) and 
    the Protocol variables."""
    def __init__(self, title, protocol, callback, master=None, hostList=['localhost'], **args):
        """ Constructor of the Form window. 
        Params:
         title: title string of the windows.
         protocol: protocol from which the form will be generated.
         callback: callback function to call when Save or Excecute are press.
        """
        Window.__init__(self, title, master, icon='scipion_bn.xbm', 
                        weight=False, minsize=(600, 450), **args)

        self.callback = callback
        self.widgetDict = {} # Store tkVars associated with params
        self.visualizeDict = args.get('visualizeDict', {})
        self.bindings = []
        self.protocol = protocol
        self.hostList = hostList
        self.protocol = protocol
        self.visualizeMode = args.get('visualizeMode', False)  # This control when to close or not after execute
        self.headerBgColor = Color.RED_COLOR
        if self.visualizeMode:
            self.headerBgColor = Color.DARK_GREY_COLOR
        self.childMode = args.get('childMode', False) # Allow to open child protocols form (for workflows)
        
        
        from pyworkflow.em import findWizards
        self.wizards = findWizards(protocol, DESKTOP_TKINTER)
        
        self.fontBig = tkFont.Font(size=12, family='helvetica', weight='bold')
        self.font = tkFont.Font(size=10, family='helvetica')#, weight='bold')
        self.fontBold = tkFont.Font(size=10, family='helvetica', weight='bold')        
        
        headerFrame = tk.Frame(self.root)
        headerFrame.grid(row=0, column=0, sticky='new')
        headerFrame.columnconfigure(0, weight=1)
        package = protocol._package
        t = '  Protocol: %s' % (protocol.getClassLabel())
        logoPath = getattr(package, '_logo', None)
        if logoPath:
            headerLabel = tk.Label(headerFrame, text=t, font=self.fontBig, image=self.getImage(logoPath), compound=tk.LEFT)
        else:
            headerLabel = tk.Label(headerFrame, text=t, font=self.fontBig)
        headerLabel.grid(row=0, column=0, padx=5, pady=(5,0), columnspan=5)
        
        def _addButton(text, icon, command, col):
            btn = tk.Label(headerFrame, text=text, image=self.getImage(icon), 
                       compound=tk.LEFT, cursor='hand2')
            btn.bind('<Button-1>', command)
            btn.grid(row=1, column=col, padx=5, sticky='se')
        
        _addButton(Message.TITLE_CITE, Icon.ACTION_REFERENCES, self._showReferences, 0)
        _addButton(Message.TITLE_DOC ,Icon.ACTION_HELP, self._showHelp, 1)
        
        if protocol.allowHeader:
            commonFrame = self._createHeaderCommons(headerFrame)
            commonFrame.grid(row=2, column=0, padx=5, pady=(0,5), 
                             sticky='news', columnspan=5)
        
        text = TaggedText(self.root, width=40, height=15, bd=0, cursor='arrow')
        text.grid(row=1, column=0, sticky='news')
        text.config(state=tk.DISABLED)
        contentFrame = self.createSections(text)
        
        bottomFrame = tk.Frame(self.root)
        bottomFrame.grid(row=2, column=0, sticky='sew')
        bottomFrame.columnconfigure(0, weight=1)
        
        btnFrame = tk.Frame(bottomFrame)
        #btnFrame.columnconfigure(2, weight=1)
        self._createButtons(btnFrame)
        btnFrame.grid(row=1, column=0, sticky='se')
        
        
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)
        
        # Resize windows to use more space if needed
        self.desiredDimensions = lambda: self.resize(contentFrame)
        #self.resize(contentFrame)
        
    def _showReferences(self, e=None):
        """ Show the list of references of the protocol. """
        self.showInfo('\n'.join(self.protocol.citations()), "References")
        
    def _showHelp(self, e=None):
        """ Show the list of references of the protocol. """
        self.showInfo(self.protocol.getDoc(), "Help")        
        
    def _createButtons(self, btnFrame):
        """ Create the bottom buttons: Close, Save and Execute. """
        btnClose = tk.Button(btnFrame, text=Message.LABEL_BUTTON_CLOSE, image=self.getImage(Icon.ACTION_CLOSE), 
                             compound=tk.LEFT, font=self.font, command=self.close)
        btnClose.grid(row=0, column=0, padx=5, pady=5, sticky='se')
        t = Message.LABEL_BUTTON_VIS
        icon = Icon.ACTION_VISUALIZE
        # Save button is not added in VISUALIZE or CHILD modes
        if not self.visualizeMode and not self.childMode:
            btnSave = tk.Button(btnFrame, text=Message.LABEL_BUTTON_RETURN, image=self.getImage(Icon.ACTION_SAVE), 
                                compound=tk.LEFT, font=self.font, command=self.save)
            btnSave.grid(row=0, column=1, padx=5, pady=5, sticky='se')
            t = Message.LABEL_BUTTON_EXEC
            icon = Icon.ACTION_EXECUTE
        # Add Execute/Visualize button
        if not self.childMode:
            btnExecute = Button(btnFrame, text=t, fg='white', bg=Color.RED_COLOR, font=self.font, 
                                image=self.getImage(icon), compound=tk.LEFT, 
                            activeforeground='white', activebackground='#A60C0C', command=self.execute)
            btnExecute.grid(row=0, column=2, padx=(5, 28), pady=5, sticky='se')
        
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
    
    def _createBoundCombo(self, parent, paramName, choices, *callbacks):
        # Create expert level combo
        param = EnumParam(choices=choices)
        var = ComboVar(param)
        #self.protocol.expertLevel.set(0)
        self._addVarBinding(paramName, var, None, *callbacks)
        #var.set(self.protocol.expertLevel.get())
        combo = ttk.Combobox(parent, textvariable=var.tkVar, 
                                state='readonly', width=10)
        combo['values'] = param.choices
        return combo
            
        
    def _createHeaderLabel(self, parent, text):
        return tk.Label(parent, text=text, font=self.font, bg='white')
    
    def _createHeaderCommons(self, parent):
        """ Create the header common values such as: runName, expertLevel, mpi... """
        commonFrame = tk.Frame(parent)
        commonFrame.columnconfigure(0, weight=1)
        commonFrame.columnconfigure(1, weight=1)
        
        ############# Create the run part ###############
        # Run name
        #runFrame = ttk.Labelframe(commonFrame, text='Run')
        runSection = SectionFrame(commonFrame, label=Message.TITLE_RUN, 
                                  headerBgColor=self.headerBgColor)
        runFrame = runSection.contentFrame
        self._createHeaderLabel(runFrame, Message.TITLE_RUN_NAME).grid(row=0, column=0, padx=5, pady=5, sticky='ne')
        entry = self._createBoundEntry(runFrame, Message.VAR_RUN_NAME, width=15, 
                                       func=self.setProtocolLabel, value=self.protocol.getObjLabel())
        entry.grid(row=0, column=1, padx=(0, 5), pady=5, sticky='nw')
        # Run Name not editable
        entry.configure(state='readonly')
        btnComment = Button(runFrame, Message.TITLE_COMMENT, Icon.ACTION_EDIT, bg='white', command=self._editObjParams)
        btnComment.grid(row=0, column=2, padx=(0, 5), pady=5, sticky='nw')
        # Run mode
        self.protocol.getDefinitionParam('')
        self._createHeaderLabel(runFrame, Message.TITLE_RUN_MODE).grid(row=1, column=0, sticky='ne', padx=5, pady=5)
        modeCombo = self._createBoundCombo(runFrame, Message.VAR_RUN_MODE, MODE_CHOICES, self._onRunModeChanged)   
        modeCombo.grid(row=1, column=1, sticky='nw', padx=(0, 5), pady=5)        
        # Expert level
        self._createHeaderLabel(runFrame, Message.TITLE_EXPERT).grid(row=2, column=0, sticky='ne', padx=5, pady=5)
        expCombo = self._createBoundCombo(runFrame, Message.VAR_EXPERT, LEVEL_CHOICES, self._onExpertLevelChanged)   
        expCombo.grid(row=2, column=1, sticky='nw', padx=(0, 5), pady=5)
        
        runSection.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        
        ############## Create the execution part ############
        # Host name
        #execFrame = ttk.Labelframe(commonFrame, text='Execution')
        execSection = SectionFrame(commonFrame, label=Message.TITLE_EXEC, 
                                  headerBgColor=self.headerBgColor)
        execFrame = execSection.contentFrame        
        self._createHeaderLabel(execFrame, Message.TITLE_EXEC_HOST).grid(row=0, column=0, padx=5, pady=5, sticky='ne')
        param = EnumParam(choices=self.hostList)
        self.hostVar = tk.StringVar()
        self._addVarBinding(Message.VAR_EXEC_HOST, self.hostVar)
        #self.hostVar.set(self.protocol.getHostName())
        expCombo = ttk.Combobox(execFrame, textvariable=self.hostVar, state='readonly', width=15)
        expCombo['values'] = param.choices        
        expCombo.grid(row=0, column=1, columnspan=3, padx=(0, 5), pady=5, sticky='nw')
        # Threads and MPI
        self._createHeaderLabel(execFrame, Message.TITLE_THREADS).grid(row=1, column=0, padx=5, pady=5, sticky='ne')
        self._createBoundEntry(execFrame, Message.VAR_THREADS).grid(row=1, column=1, padx=(0, 5))
        self._createHeaderLabel(execFrame, Message.TITLE_MPI).grid(row=1, column=2, padx=(0, 5))
        self._createBoundEntry(execFrame, Message.VAR_MPI).grid(row=1, column=3, padx=(0, 5), pady=5, sticky='nw')
        # Queue
        self._createHeaderLabel(execFrame, Message.TITLE_QUEUE).grid(row=2, column=0, padx=5, pady=5, sticky='ne', columnspan=3)
        var, frame = ParamWidget.createBoolWidget(execFrame, bg='white')
        self._addVarBinding(Message.VAR_QUEUE, var)
        frame.grid(row=2, column=3, padx=5, pady=5, sticky='nw')
        
        execSection.grid(row=0, column=1, sticky='news', padx=5, pady=5)
        
        return commonFrame
    
    def _editObjParams(self, e=None):
        """ Show a Text area to edit the protocol label and comment. """
        
        d = editObject(self, "Edit", self.root, self.protocol, self.protocol.mapper)
        
#        self.mapper = self.protocol.mapper
#        d = TextDialog(self.root, "Edit", self.protocol, self.mapper)
        
        if d.resultYes():
            label = d.valueLabel
            self.runNameVar.set(label)
            self.protocol.setObjLabel(label)
            self.protocol.setObjComment(d.valueComment)

            
    def resize(self, frame):
        self.root.update_idletasks()
        MaxHeight = 600
        MaxWidth = 600
        rh = frame.winfo_reqheight()
        rw = frame.winfo_reqwidth()
        height = min(rh + 100, MaxHeight)
        width = min(rw + 25, MaxWidth)
        x = self.root.winfo_x()
        y = self.root.winfo_y()
        self.root.geometry("%dx%d%+d%+d" % (width, height, x, y))

        return (width, height)
        
        
    def save(self, e=None):
        self._close(onlySave=True)
        
    def execute(self, e=None):
        errors = self.protocol.validate()
        
        if len(errors):
            self.showError(errors)
        else:
            self._close()
        
    def _close(self, onlySave=False):
        try:
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
        
    def _checkCondition(self, paramName):
        """Check if the condition of a param is statisfied 
        hide or show it depending on the result"""
        widget = self.widgetDict.get(paramName, None)
        
        if isinstance(widget, ParamWidget): # Special vars like MPI, threads or runName are not real widgets
            v = self.protocol.evalParamCondition(paramName) and self.protocol.evalExpertLevel(paramName)
            if v:
                widget.show()
            else:
                widget.hide()
            
    def _checkChanges(self, paramName):
        """Check the conditions of all params affected
        by this param"""
        self.setParamFromVar(paramName)
        param = self.protocol.getDefinitionParam(paramName)
        
        for d in param._dependants:
            self._checkCondition(d)
            
    def _checkAllChanges(self):
        for paramName, _ in self.protocol.iterDefinitionAttributes():
            self._checkCondition(paramName)
            
    def _onExpertLevelChanged(self, *args):
        self._checkAllChanges()
        
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
                
    def setProtocolLabel(self, paramName):
        label = self.widgetDict[paramName].get()
        self.protocol.setObjLabel(label)
        
    def setProtocolComment(self, paramName):
        label = self.widgetDict[paramName].get()
        self.protocol.setObjLabel(label)       
             
    def updateProtocolParams(self):
        for paramName, _ in self.protocol.iterDefinitionAttributes():
            self.setParamFromVar(paramName)
                
    def createSections(self, text):
        """Create section widgets"""
        r = 0
        parent = tk.Frame(text) 
        parent.columnconfigure(0, weight=1)
        tab = ttk.Notebook(parent) 
        tab.grid(row=0, column=0, sticky='news',
                 padx=5, pady=5)
        
        for section in self.protocol.iterDefinitionSections():
            label = section.getLabel()
            if label != 'General' and label != 'Parallelization':
                frame = SectionWidget(self, tab, section, 
                                      callback=self._checkChanges,
                                      headerBgColor=self.headerBgColor)
            #frame.grid(row=r, column=0, padx=10, pady=5, sticky='new')
                tab.add(frame, text=section.getLabel())
                frame.columnconfigure(0, minsize=400)
                self._fillSection(section, frame)
                r += 1
        self._checkAllChanges()
        
        # with Windows OS
        self.root.bind("<MouseWheel>", lambda e: text.scroll(e))
        # with Linux OS
        self.root.bind("<Button-4>", lambda e: text.scroll(e))
        self.root.bind("<Button-5>", lambda e: text.scroll(e)) 
        text.window_create(tk.INSERT, window=parent)
        
        return parent

def editObject(self, title, root, obj, mapper):
    """ Show a Text area to edit the protocol label and comment. """    
    return TextDialog(root, title, obj, mapper)
    

if __name__ == '__main__':
    # Just for testing
    from pyworkflow.em import ProtImportMicrographs
    p = ProtImportMicrographs()
    p.sphericalAberration.set(2.3)
    p.setSamplingRate('5.4')
    w = FormWindow(p)
    w.show()
    
   


