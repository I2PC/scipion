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
    """This class will contains the label and content frame 
    to display a Param. It will also contains a variable controlling
    the changes and the values in the GUI"""
    def __init__(self, paramName, row, label, content, var, callback=None):
        self.paramName = paramName
        self.row = row
        self.label = label
        self.content = content
        self.var = var
        self.callback = callback
        self.var.trace('w', self._onVarChanged)
        
    def _onVarChanged(self, *args):
        if self.callback is not None:
            self.callback(self.paramName)        
        
    def grid(self):
        """Grid the label and content in the specified row"""
        self.label.grid(row=self.row, column=0, sticky='ne', padx=2, pady=2)
        self.content.grid(row=self.row, column=1, padx=2, pady=2, sticky='w') 
        
    def grid_remove(self):
        self.label.grid_remove()
        self.content.grid_remove()
        
    def set(self, value):
        self.var.set(value)
        
    def get(self):
        return self.var.get()
        

    
class FormWindow(Window):
    """This class will create the Protocol params GUI to enter parameters.
    The creaation of compoments will be based on the Protocol Form definition.
    This class will serve as a connection between the GUI variables (tk vars) and 
    the Protocol variables."""
    def __init__(self, title, protocol, master=None, **args):
        Window.__init__(self, title, master, weight=False, **args)

        self.varDict = {} # Store tkVars associated with params
        
        self.fontBig = tkFont.Font(size=12, family='verdana', weight='bold')
        self.font = tkFont.Font(size=10, family='verdana')#, weight='bold')
        self.fontBold = tkFont.Font(size=10, family='verdana', weight='bold')
        
        headerFrame = tk.Frame(self.root)
        headerFrame.grid(row=0, column=0, sticky='new')
        headerLabel = tk.Label(headerFrame, text='Protocol: Import micrographs', font=self.fontBig)
        headerLabel.grid(row=0, column=0, padx=5, pady=5)
        
        text = TaggedText(self.root, width=40, height=15, bd=0, cursor='arrow')
        text.grid(row=1, column=0, sticky='news')
        text.config(state=tk.DISABLED)
        self.protocol = protocol
        self.createSections(text)
        
        btnFrame = tk.Frame(self.root)
        btnFrame.columnconfigure(0, weight=1)
        btnFrame.grid(row=2, column=0, sticky='sew')
        # Add create project button
        btn = tk.Button(btnFrame, text='Execute', fg='white', bg='#7D0709', font=self.font, 
                        activeforeground='white', activebackground='#A60C0C', command=self.execute)
        btn.grid(row=0, column=0, padx=(5, 100), pady=5, sticky='se')
        
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)
        
    def execute(self, e=None):
        self.updateProtocolParams()
        self.close()
        self.protocol.run()
    
           
    def _fillSection(self, sectionParam, sectionFrame):
        parent = sectionFrame.contentFrame
        r = 0
        for paramName, param in sectionParam.iterParams():
            if sectionParam.questionParam.get() != paramName: # Skip if questionParam was set in section
                protVar = getattr(self.protocol, paramName, None)
                if protVar is None:
                    raise Exception("_fillSection: param '%s' not found in protocol" % paramName)
                # Create the label
                if param.isImportant:
                    f = self.fontBold
                else:
                    f = self.font
                label = tk.Label(parent, text=param.label.get(), bg='white', font=f)
                # Create the content to display and the variable associated
                content, var = self._createWidgetFromParam(parent, param)
                # Set the value to the variable from Protocol instance or default
                var.set(protVar.get(param.default.get()))
                self.varDict[paramName] = ParamWidget(paramName, r, label, content, var, self._checkChanges)
                self._checkCondition(paramName)
                r += 1
            
    def _createWidgetFromParam(self, parent, param):
        """Create the widget for an specific param.
        The function should return a tuple (content, var)
        content: will be the outer container of the widget(tipically a tk.Frame
        var: should contain .set and .get methods with the same type of the protocol var"""
        # Create widgets for each type of param
        t = type(param)
        #TODO: Move this to a Renderer class to be more flexible
        if t is BooleanParam:
            content = tk.Frame(parent, bg='white')
            var = BoolVar()
            rb1 = tk.Radiobutton(content, text='Yes', bg='white', variable=var.tkVar, value=1)
            rb1.grid(row=0, column=0, padx=2)
            rb2 = tk.Radiobutton(content, text='No', bg='white', variable=var.tkVar, value=0)
            rb2.grid(row=0, column=1, padx=2)
            
        elif t is EnumParam:
            var = tk.IntVar()
            if param.display == EnumParam.DISPLAY_COMBO:
                combo = ttk.Combobox(parent, textvariable=var, state='readonly')
                combo['values'] = param.choices
                content = combo
            elif param.display == EnumParam.DISPLAY_LIST:
                rbFrame = tk.Frame(parent)
                for i, opt in enumerate(param.choices):
                    rb = tk.Radiobutton(rbFrame, text=opt, variable=var, value=i)
                    rb.grid(row=i, column=0, sticky='w')
                content = rbFrame
            else:
                raise Exception("Invalid display value '%s' for EnumParam" % str(param.display))
        elif t is PointerParam:
            var = PointerVar()
            content = tk.Entry(parent, width=25, textvariable=var.tkVar)
        else:
            #v = self.setVarValue(paramName)
            var = tk.StringVar()
            content = tk.Entry(parent, width=25, textvariable=var)
                
            
        return (content, var)
        
    def _checkCondition(self, paramName):
        """Check if the condition of a param is statisfied 
        hide or show it depending on the result"""
        w = self.varDict[paramName]
        print " checking condition for: ", paramName
        v = self.protocol._definition.evalCondition(self.protocol, paramName)
        print " result: ", v
        if v:
            w.grid()
        else:
            w.grid_remove()
            
    def _checkChanges(self, paramName):
        """Check the conditions of all params affected
        by this param"""
        self.setParamFromVar(paramName)
        param = self.protocol._definition.getParam(paramName)
        
        #TODO: REMOVE DEBUG PRINTS
        #print "_checkChanges"
        #print "param changed: ", paramName
        #print " dependants: ", param._dependants
        
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
            w = self.varDict[paramName]
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
        
  
if __name__ == '__main__':
    # Just for testing
    from pyworkflow.em import ProtImportMicrographs
    p = ProtImportMicrographs()
    p.sphericalAberration.set(2.3)
    p.samplingRate.set('5.4')
    w = FormWindow(p)
    w.show()
    
   


