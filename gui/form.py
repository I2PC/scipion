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
from protocol.params import PointerParam
"""
This modules implements the automatic
creation of protocol form GUI from its
params definition.
"""

import Tkinter as tk
import ttk
import tkFont

import gui
from gui import configureWeigths, Window
from text import TaggedText
from widgets import Button
from dialog import showInfo
from pyworkflow.protocol.params import *


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
        if value:
            v = '%s.%s' % (value.getName(), value.strId())
        self.tkVar.set(v)   
            
    def get(self):
        return self.value        
        
class SectionFrame(tk.Frame):
    """This class will be used to create a section in FormWindow"""
    def __init__(self, form, master, section, **args):
        tk.Frame.__init__(self, master, **args)
        self.form = form
        self.section = section
        self.__createHeader()
        self.__createContent()
        
    def __createHeader(self):
        bgColor = gui.cfgButtonBgColor
        self.headerFrame = tk.Frame(self, bd=2, relief=tk.RAISED, bg=bgColor)
        self.headerFrame.grid(row=0, column=0, sticky='new')
        configureWeigths(self.headerFrame)
        self.headerFrame.columnconfigure(1, weight=1)
        #self.headerFrame.columnconfigure(2, weight=1)
        self.headerLabel = tk.Label(self.headerFrame, text=self.section.label.get(), fg='white', bg=bgColor)
        self.headerLabel.grid(row=0, column=0, sticky='nw')
        
        if self.section.hasQuestion():
            question = self.section.getQuestion()             
            self.tkVar = BoolVar()
            self.tkVar.set(question.get())
            
            self.chbLabel = tk.Label(self.headerFrame, text=question.label.get(), fg='white', bg=bgColor)
            self.chbLabel.grid(row=0, column=1, sticky='e', padx=2)
            
            self.chb = tk.Checkbutton(self.headerFrame, variable=self.tkVar.tkVar, 
                                      bg=bgColor, activebackground=gui.cfgButtonActiveBgColor, #bd=0,
                                      command=self.__expandCollapse)
            self.chb.grid(row=0, column=2, sticky='e')
                    #bg=SectionBgColor, activebackground=ButtonActiveBgColor)        
    
    def __createContent(self):
        self.contentFrame = tk.Frame(self, bg='white', bd=0)
        self.contentFrame.grid(row=1, column=0, sticky='news', padx=5, pady=5)
        configureWeigths(self.contentFrame)
        self.columnconfigure(0, weight=1)
        
    def __expandCollapse(self, e=None):
        if self.tkVar.get():
            self.contentFrame.grid(row=1, column=0, sticky='news', padx=5, pady=5)
        else:
            self.contentFrame.grid_remove()
            
    def get(self):
        """Return boolean value if is selected"""
        return self.tkVar.get() == 1
    
    def set(self, value):
        pass
    
    
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
    def __init__(self, row, paramName, param, window, parent, value, callback=None):
        self.window = window
        self.row = row
        self.paramName = paramName
        self.param = param
        self.parent = parent
        
        self._createLabel() # self.label should be set after this 
        self._btnCol = 0
        self._createButtonsFrame() # self.btnFrame should be set after this
        self._createContent() # self.content and self.var should be set after this
        
        self.set(value)
        self.callback = callback
        self.var.trace('w', self._onVarChanged)
        
        
    def _createLabel(self):
        if self.param.isImportant.get():
            f = fontBold
        else:
            f = font
        self.label = tk.Label(self.parent, text=self.param.label.get(), 
                              bg='white', font=f)
               
    def _createContent(self):
        self.content = tk.Frame(self.parent, bg='white')
        gui.configureWeigths(self.content)
        self._createContentWidgets(self.param, self.content) # self.var should be set after this
        
    def _addButton(self, text, imgPath, cmd):
        btn = Button(self.btnFrame, text, imgPath, bg='white', command=cmd)
        btn.grid(row=0, column=self._btnCol, sticky='w', padx=2)
        self._btnCol += 1
        
    def _createButtonsFrame(self):
        self.btnFrame = tk.Frame(self.parent, bg='white')
        if self.param.help.hasValue():
            self._addButton('Help', 'contents24.png', self._showHelpMessage)
        
    def _showHelpMessage(self, e=None):
        showInfo("Help", self.param.help.get(), self.parent)
        
    def _createContentWidgets(self, param, content):
        """Create the specific widgets inside the content frame"""
        # Create widgets for each type of param
        t = type(param)
        #TODO: Move this to a Renderer class to be more flexible
        if t is BooleanParam:
            var = BoolVar()
            frame = tk.Frame(content, bg='white')
            frame.grid(row=0, column=0, sticky='w')
            rb1 = tk.Radiobutton(frame, text='Yes', bg='white', variable=var.tkVar, value=1)
            rb1.grid(row=0, column=0, padx=2, sticky='w')
            rb2 = tk.Radiobutton(frame, text='No', bg='white', variable=var.tkVar, value=0)
            rb2.grid(row=0, column=1, padx=2, sticky='w')
            
        elif t is EnumParam:
            var = tk.IntVar()
            if param.display == EnumParam.DISPLAY_COMBO:
                combo = ttk.Combobox(content, textvariable=var, state='readonly')
                combo['values'] = param.choices
                combo.grid(row=0, column=0)
            elif param.display == EnumParam.DISPLAY_LIST:
                for i, opt in enumerate(param.choices):
                    rb = tk.Radiobutton(content, text=opt, variable=var, value=i)
                    rb.grid(row=i, column=0, sticky='w')
            else:
                raise Exception("Invalid display value '%s' for EnumParam" % str(param.display))
        
        elif t is PointerParam:
            var = PointerVar()
            entry = tk.Entry(content, width=25, textvariable=var.tkVar)
            entry.grid(row=0, column=0)
            self._addButton("Select", 'zoom.png', self._browseObject)
        else:
            #v = self.setVarValue(paramName)
            var = tk.StringVar()
            entryWidth = 25
            if t is FloatParam or t is IntParam:
                entryWidth = 10 # Reduce the entry width for numbers entries
            entry = tk.Entry(content, width=entryWidth, textvariable=var)
            entry.grid(row=0, column=0, sticky='w')
            
        self.var = var
        
    def _browseObject(self, e=None):
        """Select an object from DB
        This function is suppose to be used only for PointerParam"""
        from pyworkflow.gui.dialog import SubclassesTreeProvider, ListDialog
        tp = SubclassesTreeProvider(self.window.protocol.mapper, self.param.pointerClass.get())
        dlg = ListDialog(self.parent, "Select object", tp)
        if dlg.value is not None:
            self.set(dlg.value)
        
    def _onVarChanged(self, *args):
        if self.callback is not None:
            self.callback(self.paramName)        
        
    def grid(self):
        """Grid the label and content in the specified row"""
        self.label.grid(row=self.row, column=0, sticky='ne', padx=2, pady=2)
        self.content.grid(row=self.row, column=1, padx=2, pady=2, sticky='news')
        self.btnFrame.grid(row=self.row, column=2, padx=2, sticky='nw')
        
    def grid_remove(self):
        self.label.grid_remove()
        self.content.grid_remove()
        self.btnFrame.grid_remove()
        
    def set(self, value):
        if value is not None:
            self.var.set(value)
        
    def get(self):
        return self.var.get()
        

# FIXME: Organize this and make it configurable
fontBig = None
font = None
fontBold = None
    
class FormWindow(Window):
    """This class will create the Protocol params GUI to enter parameters.
    The creaation of compoments will be based on the Protocol Form definition.
    This class will serve as a connection between the GUI variables (tk vars) and 
    the Protocol variables."""
    def __init__(self, title, protocol, executeCallback, master=None, **args):
        Window.__init__(self, title, master, icon='scipion_bn.xbm', weight=False, **args)

        self.executeCallback = executeCallback
        self.widgetDict = {} # Store tkVars associated with params
        
        global fontBig, font, fontBold
        fontBig = tkFont.Font(size=12, family='verdana', weight='bold')
        font = tkFont.Font(size=10, family='verdana')#, weight='bold')
        fontBold = tkFont.Font(size=10, family='verdana', weight='bold')        

        
        headerFrame = tk.Frame(self.root)
        headerFrame.grid(row=0, column=0, sticky='new')
        headerLabel = tk.Label(headerFrame, text='Protocol: ' + protocol.getClassName(), font=fontBig)
        headerLabel.grid(row=0, column=0, padx=5, pady=5)
        
        text = TaggedText(self.root, width=40, height=15, bd=0, cursor='arrow')
        text.grid(row=1, column=0, sticky='news')
        text.config(state=tk.DISABLED)
        self.protocol = protocol
        contentFrame = self.createSections(text)
        
        bottomFrame = tk.Frame(self.root)
        bottomFrame.columnconfigure(0, weight=1)
        bottomFrame.grid(row=2, column=0, sticky='sew')
        
        expertFrame = tk.Frame(bottomFrame)
        #expertFrame.columnconfigure(0, weight=1)
        expertFrame.grid(row=0, column=0, sticky='sw', padx=5, pady=5)
        expLabel = tk.Label(expertFrame, text="Expert level: ", font=font)
        expLabel.grid(row=0, column=0)
        var = tk.StringVar()
        var.set(LEVEL_CHOICES[0])
        expCombo = ttk.Combobox(expertFrame, textvariable=var, state='readonly')
        expCombo['values'] = LEVEL_CHOICES
        expCombo.grid(row=0, column=1)
        self.expertVar = var
        
        btnFrame = tk.Frame(bottomFrame)
        btnFrame.columnconfigure(2, weight=1)
        btnFrame.grid(row=1, column=0, sticky='sew')
        # Save and close buttons
        def showMsg(msg):
            print "clicked: ", msg
        
        btnClose = tk.Button(btnFrame, text="Close", image=self.getImage('dialog_close.png'), compound=tk.LEFT, font=font,
                          command=self.close)
        btnClose.grid(row=0, column=0, padx=5, pady=5, sticky='sw')
        btnSave = tk.Button(btnFrame, text="Save", image=self.getImage('filesave.png'), compound=tk.LEFT, font=font, 
                          command=lambda : showMsg('Save'))
        btnSave.grid(row=0, column=1, padx=5, pady=5, sticky='sw')
        # Add create project button
        
        btnExecute = Button(btnFrame, text='   Execute   ', fg='white', bg='#7D0709', font=font, 
                        activeforeground='white', activebackground='#A60C0C', command=self.execute)
        btnExecute.grid(row=0, column=2, padx=(15, 50), pady=5, sticky='se')
        
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)
        
        # Resize windows to use more space if needed
        self.desiredDimensions = lambda: self.resize(contentFrame)
        #self.resize(contentFrame)
        
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
        
    def execute(self, e=None):
        self.close()
        self.executeCallback(self.protocol)
           
    def _fillSection(self, sectionParam, sectionFrame):
        parent = sectionFrame.contentFrame
        r = 0
        for paramName, param in sectionParam.iterParams():
            if sectionParam.questionParam.get() != paramName: # Skip if questionParam was set in section
                protVar = getattr(self.protocol, paramName, None)
                if protVar is None:
                    raise Exception("_fillSection: param '%s' not found in protocol" % paramName)
                # Create the label
                self.widgetDict[paramName] = ParamWidget(r, paramName, param, self, parent, 
                                                         value=protVar.get(param.default.get()),
                                                         callback=self._checkChanges)
                # Set the value to the variable from Protocol instance or default
                self._checkCondition(paramName)
                r += 1         
            
        
    def _checkCondition(self, paramName):
        """Check if the condition of a param is statisfied 
        hide or show it depending on the result"""
        widget = self.widgetDict[paramName]
        v = self.protocol._definition.evalCondition(self.protocol, paramName)
        if v:
            widget.grid()
        else:
            widget.grid_remove()
            
    def _checkChanges(self, paramName):
        """Check the conditions of all params affected
        by this param"""
        self.setParamFromVar(paramName)
        param = self.protocol._definition.getParam(paramName)
        
        for d in param._dependants:
            self._checkCondition(d)
        
    def getVarValue(self, varName):
        """This method should retrieve a value from """
        pass
    
    def setVarFromParam(self, tkVar, paramName):
        value = getattr(self.protocol, paramName, None)
        if value is not None:
            tkVar.set(value.get(''))
           
    def setParamFromVar(self, paramName):
        value = getattr(self.protocol, paramName, None)
        if value is not None:
            w = self.widgetDict[paramName]
            value.set(w.get())
             
    def updateProtocolParams(self):
        for paramName, _ in self.protocol._definition.iterParams():
            self.setParamFromVar(paramName)
                
    def createSections(self, text):
        """Load the list of projects"""
        r = 0
        parent = tk.Frame(text) 
        parent.columnconfigure(0, weight=1)
        
        for section in self.protocol._definition.iterSections():
            frame = SectionFrame(self, parent, section, bg='white')
            frame.grid(row=r, column=0, padx=10, pady=5, sticky='new')
            frame.columnconfigure(0, minsize=300)
            self._fillSection(section, frame)
            r += 1
        
        # with Windows OS
        self.root.bind("<MouseWheel>", lambda e: text.scroll(e))
        # with Linux OS
        self.root.bind("<Button-4>", lambda e: text.scroll(e))
        self.root.bind("<Button-5>", lambda e: text.scroll(e)) 
        text.window_create(tk.INSERT, window=parent)
        
        return parent
  
if __name__ == '__main__':
    # Just for testing
    from pyworkflow.em import ProtImportMicrographs
    p = ProtImportMicrographs()
    p.sphericalAberration.set(2.3)
    p.samplingRate.set('5.4')
    w = FormWindow(p)
    w.show()
    
   


