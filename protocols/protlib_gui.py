#!/usr/bin/env xmipp_python
'''
#/***************************************************************************
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
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
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************
 '''

import os
from glob import glob
from subprocess import Popen
from os.path import join, relpath, exists
import Tkinter as tk
import tkFont


from protlib_base import getWorkingDirFromRunName, getExtendedRunName
from protlib_utils import loadModule, which, runShowJ
from protlib_gui_ext import centerWindows, changeFontSize, askYesNo, Fonts, registerCommonFonts, \
    showError, showInfo, showBrowseDialog, showWarning, AutoScrollbar, FlashMessage,\
    TaggedText, XmippText, XmippButton, OutputText
from protlib_filesystem import getXmippPath, xmippExists, removeBasenameExt
from config_protocols import *
from protlib_sql import SqliteDb
from protlib_include import *
from protlib_parser import ProtocolParser
from protlib_xmipp import redStr, greenStr


class ProtocolStyle():
    ''' Class to define some style settings like font, colors, etc '''
    def __init__(self, configModuleName=None):        
        #Font
        self.FontName = FontName
        self.FontSize = FontSize
        self.ButtonFontSize = self.FontSize
        #TextColor
        self.CitationTextColor = "dark olive green"
        self.LabelTextColor = "black"
        self.SectionTextColor = "blue4"
        #Background Color
        self.BgColor = "white"
        self.LabelBgColor = self.BgColor
        self.HighlightBgColor = self.BgColor
        self.ButtonBgColor = "LightBlue"
        self.ButtonActiveBgColor = "LightSkyBlue"
        self.EntryBgColor = "lemon chiffon" 
        self.ExpertLabelBgColor = "light salmon"
        #Color
        self.ListSelectColor = "DeepSkyBlue4"
        self.BooleanSelectColor = "DeepSkyBlue4"
        #Dimensions limits
        self.MaxHeight = 600
        self.MaxWidth = 800
        self.MaxFontSize = 14
        self.MinFontSize = 6
        self.WrapLenght = self.MaxWidth / 2
        
        if configModuleName:
            self.load(configModuleName)
                    
    def load(self, configModuleName):
        mod = loadModule(configModuleName, False)
        if mod:
            modDir = dir(mod)
            selfDir = dir(self)
            for a in modDir:
                if a in selfDir and not a.startswith('_'):
                        self.__dict__[a] = mod.__dict__[a]
                    
    def createFonts(self):
        self.Font = tkFont.Font(family=self.FontName, size=self.FontSize, weight=tkFont.BOLD)

def createSection(parent, text):
    frame = tk.Frame(parent, bd=2, relief=tk.RAISED, bg=SectionBgColor)
    frame.columnconfigure(0, weight=1)
    label = tk.Label(frame, text=text, fg=LabelTextColor, bg=SectionBgColor, font=Fonts['button'])
    label.grid(row=0, column=0, sticky=tk.W)        
    content = tk.Frame(frame, bg=LabelBgColor, bd=0)
    content.grid(row=1, column=0, columnspan=5, sticky=tk.NSEW, ipadx=5, ipady=5)  
    return (frame, label, content)     

class BasicGUI(): 
    
    def createBasicGUI(self, master=None):
        """Perform some basic GUI initializations.
        - create the main window
        - create style and fonts
        - create scrollable canvas
        """
        if not master:
            master = tk.Tk()
        self.master = master
        self.master.withdraw()
        self.style = ProtocolStyle('config_protocols')
        self.style.createFonts()
        
    def createScrollableCanvas(self):
        """Create an scrollable canvas
        It will set .frame and .canvas new properties
        """
        vscrollbar = AutoScrollbar(self.master)
        vscrollbar.grid(row=0, column=1, sticky='ns')
        hscrollbar = AutoScrollbar(self.master, orient=tk.HORIZONTAL)
        hscrollbar.grid(row=1, column=0, sticky='ew')
        self.canvas = tk.Canvas(self.master, background=self.style.BgColor,
                        yscrollcommand=vscrollbar.set,
                        xscrollcommand=hscrollbar.set)
        self.canvas.grid(row=0, column=0, sticky='nsew')
        vscrollbar.config(command=self.canvas.yview)
        hscrollbar.config(command=self.canvas.xview)
        self.master.grid_rowconfigure(0, weight=1)
        self.master.grid_columnconfigure(0, weight=1)
        self.frame = tk.Frame(self.canvas, background=self.style.BgColor)
        self.frame.rowconfigure(0, weight=1)
        self.frame.columnconfigure(0, weight=1)
        
    def resize(self):
        height = min(self.frame.winfo_reqheight() + 25, MaxHeight)
        width = min(self.frame.winfo_reqwidth() + 25, MaxWidth)
        x = self.master.winfo_x()
        y = self.master.winfo_y()
        self.master.geometry("%dx%d%+d%+d" % (width, height, x, y))
        return (width, height)

    def updateScrollRegion(self):
        self.frame.update_idletasks()
        self.canvas.config(scrollregion=self.canvas.bbox("all")) 
        
    def launchCanvas(self):
        # Launch the window
        self.canvas.create_window(0, 0, anchor='nw', window=self.frame)
        self.updateScrollRegion()
    
    def launchGUI(self):
        self.launchCanvas() 
        centerWindows(self.master, self.resize() )
        self.master.deiconify()
        self.master.lift()
        self.master.mainloop() 

    def currentRow(self, parent=None):
        if not parent:
            parent = self.frame
        return parent.grid_size()[1]
    
    def nextRow(self, parent=None):
        if not parent:
            parent = self.frame
        return parent.grid_size()[1] + 1
    
    def addLabel(self, text, row=None, column=0, columnspan=1, sticky='ew', bgColor='', fgColor='', parent=None):
        if not parent:
            parent = self.frame
        if fgColor == '':
            fgColor = self.style.LabelTextColor
        if bgColor == '':
            bgColor = self.style.LabelBgColor
        if not row:
            row = self.nextRow(parent)
        label = tk.Label(parent, text=text, bg=bgColor, fg=fgColor)
        label.grid(row=row, column=column, sticky=sticky, columnspan=columnspan)
        return label
    
    def addCheckButton(self, text, row, column, controlVar, default, command, sticky, parent=None):
        if not parent:
            parent = self.frame
        controlVar.set(default)
        check = tk.Checkbutton(parent, text=text, variable=controlVar,
                      command=command,
                      selectcolor=self.style.BooleanSelectColor,
                      bg=self.style.BgColor,
                      font=Fonts['normal'])
        check.grid(row=row, column=column, sticky=sticky)
        return check

    def addRadioButton(self, text, row, column, variable, value, command, sticky='', parent=None):
        if not parent:
            parent = self.frame
        radio = tk.Radiobutton(parent, text=text, variable=variable,
                          value=value, indicatoron=0,
                          command=command,
                          bg=self.style.ButtonBgColor,
                          activebackground=self.style.ButtonActiveBgColor,
                          highlightbackground=self.style.HighlightBgColor,
                          selectcolor=self.style.ButtonActiveBgColor,
                          font=Fonts['normal'])
        radio.grid(row=row, column=column, sticky=sticky)
        return radio

    def addButton(self, text, row, column, command, sticky="", binding="", parent=None):
        if not parent:
            parent = self.frame
        button = tk.Button(parent, text=text,
                        command=command,
                        bg=self.style.ButtonBgColor,
                        activebackground=self.style.ButtonActiveBgColor)
        button.grid(row=row, column=column, sticky=sticky)
        if binding != "":
            self.master.bind(binding, command)
        return button

    def addLine(self, _color, column, columnspan, parent=None):
        if not parent:
            parent = self.frame
        line = tk.Frame(parent, height=2, bd=1, bg=_color, relief=tk.RIDGE)
        line.grid(row=self.nextRow(), column=column, columnspan=columnspan, sticky='ew')
        return line
        

class ProtocolWidget():
    def __init__(self, master, var):
        self.master = master
        self.widgetslist = [] # Store real tk widgets
        self.childwidgets = [] # This will only be used for sections
        self.var = var
       
    def satisfiesCondition(self):
        if self.var:
            return self.var.satisfiesCondition()
        else:
            return True
     
    def checkVisibility(self):
        '''Show or hide the widget, also return the visibility'''
        show = self.getVisibility()
        self.display(show)
        return show    
        
    def getVisibility(self):
        if self.var.isHidden():
            return False
        show = self.var.isVisualize() == self.master.visualize_mode
        if show and self.var.isExpert():
            show = self.master.expert_mode
        if show:
            show = self.satisfiesCondition() #has condition, check it
        return show
        
            
    def display(self, value):
        for w in self.widgetslist:
            if value:
                w.grid()
            else:
                w.grid_remove()
                
    def validate(self):
        errors = []
        wiz = loadModule("protlib_wizard")
        if self.getVisibility():
            for v in self.var.validators:
                if hasattr(wiz, v):
                    e = getattr(wiz, v)(self.var)
                else:
                    e = "Validator <%s> for <%s> not found" % (v, self.var.comment)
                if e:
                    errors.append(e)
        return errors
    
    def isSectionExpanded(self):
        ''' Return true if the section is visible and expanded'''
        if not self.var.hasQuestion():
            return True
        return self.tkvar.get() == 'True' and self.getVisibility()

class ProtocolGUI(BasicGUI):
        
    def initialize(self, project, run, saveCallback):
        ''' Internal function to initialize some properties.
        It is not the constructor '''
        self.expert_mode = True
        self.visualize_mode = False
        self.lastrow = 0
        self.sectionslist = [] # List of all sections
        self.parser = None
        self.run = run
        self.project = project
        # Script title
        self.programname = removeBasenameExt(self.run['source']) 
        self.maxLabelWidth = 0
        self.hasVisualizeOptions = False
        #self.master = None
        self.Error = None
        self.saveCallback = saveCallback
        registerCommonFonts()

    #-------------------------------------------------------------------
    # Widgets creation and GUI building
    #-------------------------------------------------------------------  
    def addButton(self, text, cmd, underline, row, col, sticky, imageFilename=None, parent=None, tip=None):
        f = Fonts['button']
        helpImage = None
        if parent is None:
            parent = self.frame
        if imageFilename:
            try:
                imgPath = join(getXmippPath('resources'), imageFilename)
                helpImage = tk.PhotoImage(file=imgPath)
            except tk.TclError:
                pass
        
        if helpImage:
            btn = tk.Button(parent, image=helpImage, bg=self.style.LabelBgColor, bd=0)
            btn.image = helpImage
            pad = 3
        else:
            btn = tk.Button(parent, text=text, underline=underline, font=f,
                     bg=self.style.ButtonBgColor)
            pad = 5
        btn.config(command=cmd, activebackground=self.style.ButtonActiveBgColor)
        btn.grid(row=row, column=col, sticky=sticky, padx=pad, pady=pad)
        if tip:
            from protlib_gui_ext import ToolTip
            ToolTip(btn, tip, 500)
        return btn
    
    def addRadioButton(self, w, var, text, value, row, col, parent):
        rb = tk.Radiobutton(parent, text=text, variable=var.tkvar, value=value, bg=self.style.LabelBgColor,
                                 font=Fonts['normal'])
        rb.grid(row=row, column=col, sticky='w')
        w.widgetslist.append(rb) 
        return rb
        
    def expandCollapseSection(self, section):
        '''Expand/collapse the section area, return True if expanded'''
        if section.tkvar.get() == 'False':
            section.content.grid_remove()
            expanded = False
        else:
            for w in section.childwidgets:
                w.checkVisibility()
            section.content.grid(row=1, column=0, columnspan=5, sticky='nsew')
            expanded = True
        self.updateScrollRegion()
        return expanded

    def createSectionWidget(self, var):
        w = ProtocolWidget(self, var)
        w.frame, w.label, w.content = createSection(self.frame, var.comment)
        w.name = var.comment
        w.widgetslist.append(w.frame)
        w.frame.grid(row=self.getRow(), column=0, columnspan=5, 
                         sticky='ew', pady=5, padx=(10,5))
        self.sectionslist.append(w)
        return w
    
    def createWidget(self, section, var):
        w = ProtocolWidget(self, var)  
        
        if var.isHidden():
            return w
        
        label_row = row = self.getRow() # Get row where to place the widget        
        label_text = var.comment
        label_color = self.style.LabelTextColor
        label_bgcolor = self.style.LabelBgColor        
        
        section.childwidgets.append(w)
        frame = section.content
        #widgets inherit expert from section and its conditions 
        if section.var.isExpert():
            var.setTag('expert')
        if section.var.isVisualize():
            self.hasVisualizeOptions = True
            var.setTag('visualize')
        
        var.updateValidators()
        
        if var.isExpert():
            label_bgcolor = self.style.ExpertLabelBgColor
            
        if var.isText():
            var.tktext = tk.Text(frame, width=66, height=10, wrap=tk.WORD, bg=EntryBgColor, font=Fonts['normal'])
            var.tktext.grid(row=label_row, column=0, columnspan=5, sticky='ew', padx=(10, 0), pady=(10, 0))
            var.tktext.insert(tk.END, var.value)
            w.widgetslist.append(var.tktext)
            return w

        var_column = 1
        var.tkvar = tk.StringVar()   
        var.setTkValue(var.getValue())
        keys = var.tags.keys()
            
        if var.isBoolean():
            #Check if that boolean variable is the question of the section
            if section.var.hasQuestion() and len(section.childwidgets) == 1:
                section.tkvar = var.tkvar
                #Label(section.frame, text=label_text, bg=SectionBgColor).grid(row=0, column=1, padx=(5, 0))
                chb = tk.Checkbutton(section.frame, text=label_text, variable=var.tkvar, 
                            onvalue='True', offvalue='False', font=Fonts['normal'],
                            command=lambda:self.expandCollapseSection(section),
                            bg=SectionBgColor, activebackground=ButtonActiveBgColor)
                chb.grid(row=0, column=1, padx=(5, 0))
                self.expandCollapseSection(section)
                return w
            else:
                self.addRadioButton(w, var, 'Yes', 'True', row, var_column, frame)  
                self.addRadioButton(w, var, 'No', 'False', row, var_column + 1, frame) 
                if 'view' in keys:
                    btn = self.addButton("View", lambda:self.visualizeVar(var.name), -1, row, var_column+2, 'nw', 'visualize.gif', frame, 'View')
                    w.widgetslist.append(btn)
                var.tkvar.trace('w', self.checkVisibility)
        elif 'list_combo' in keys:
            items = var.getTagValues('list_combo')
            optMenu = tk.OptionMenu(frame, var.tkvar, *items)
            optMenu.config(bg=ButtonBgColor, activebackground=ButtonActiveBgColor, font=Fonts['normal'])
            optMenu['menu'].config(bg=ButtonBgColor, activebackground=ButtonActiveBgColor, font=Fonts['normal'])
            optMenu.grid(row=row, column=var_column, sticky='ew', columnspan=2)
            w.widgetslist.append(optMenu)
            var.tkvar.trace('w', self.checkVisibility)            
        elif 'list' in keys:
            items = var.tags['list'].split(',')
            for o in items:
                o = o.strip()
                self.addRadioButton(w, var, o, o, row, var_column, frame)
                row = self.getRow()
            var.tkvar.trace('w', self.checkVisibility)
         
        else: #Add a text Entry
            from protlib_gui_ext import AutoCompleteEntry
            entry = AutoCompleteEntry(frame, textvariable=var.tkvar, bg=EntryBgColor, font=Fonts['normal'])
            entry.grid(row=row, column=var_column, columnspan=2, sticky='ew')
            w.widgetslist.append(entry)
            args = None
            
            def getEntries(var, onlyDir=False):
                cwd = os.getcwd()
                pattern = join(cwd, var.tkvar.get()) + '*'
                entries = []
                for p in glob(pattern):
                    p = relpath(p, cwd)
                    if os.path.isdir(p):
                        p+="/"
                        entries.append(p)
                    else:
                        if not onlyDir:
                            entries.append(p)
                return entries
            
            wiz = loadModule('protlib_wizard')
            viewFunc = None

            def getRuns(var):
                return self.project.getFinishedRunList(var.getTagValues('run'))

            def getWizardFunc(wizName, default=wiz.wizardNotFound):
                return getattr(wiz, wizName, default)
                
            if 'file' in keys:
                viewFunc = wiz.wizardShowJ
                entry.setBuildListFunction(lambda: getEntries(var), refreshOnTab=True)
                args = ['Browse', lambda: wiz.wizardBrowse(self, var), 'fileopen.gif', 'Browse file']
            elif 'dir' in keys:
                entry.setBuildListFunction(lambda: getEntries(var,onlyDir=True), refreshOnTab=True)
                args = ['Browse', lambda: wiz.wizardBrowse(self, var), 'folderopen.gif', 'Browse folder']
            elif 'run' in keys:
                entry.setBuildListFunction(lambda: getRuns(var))
                runList = getRuns(var)
                func = None
                if 'post_run' in keys:
                    func = getWizardFunc(var.tags['post_run'], None)
                if len(runList) == 1:
                    var.setTkValue(runList[0])
                    if func is not None:
                        func(self)
                args = ['Select Run', lambda: self.selectFromList(var, runList, func), 'wizard.gif', 'Select run']
                # Run are always input variables and should not be empty
                var.validators.append('validatorNonEmpty')
            elif 'blocks' in keys:
                #md = self.varsDict[var.tags['blocks']].getValue()
                args = ['Select Blocks', lambda: self.selectFromList(var, ['block1', 'block2', 'block3']), 'wizard.gif', 'Select blocks']
            elif 'wizard' in keys:
                funcName = var.tags['wizard']
                func = getWizardFunc(funcName)
                args = [funcName, lambda:func(self, var), 'wizard.gif', funcName]
            if args:
                btn = self.addButton(args[0], args[1], -1, label_row, var_column+2, 'nw', args[2], frame, args[3])
                w.widgetslist.append(btn)
            if 'view' in keys:
                viewFunc = wiz.wizardShowJ
                funcName = var.tags['view']
                viewFunc = getattr(wiz, funcName, viewFunc)
            if viewFunc:
                btn = self.addButton("View", lambda:viewFunc(self, var), -1, label_row, var_column+3, 'nw', 'visualize.gif', frame, 'View')
                w.widgetslist.append(btn)
        if var.help:
            btn = self.addButton("Help", lambda: self.showHelp(var.help.replace('"', '')), -1, label_row, var_column+4, 'nw', 'help.gif', frame, 'Show info')
            w.widgetslist.append(btn)
        if var.name == 'RunName':
            label_text += ' %s_' % self.run['protocol_name']
        #label_text="_%s_"%label_text
        label = tk.Label(frame, text=label_text, fg=label_color, bg=label_bgcolor, font=Fonts['normal'])
        label.grid(row=label_row, column=0, sticky='e', padx=(5, 10))
        
        self.maxLabelWidth = max(self.maxLabelWidth, label.winfo_reqwidth())
        w.widgetslist.append(label)
    
        return w
        
    def visualizeVar(self, varName):
        prot = self.getProtocol()        
        prot.visualizeVar(varName)
        
    def setTitleAndHeader(self):
        self.titleText = "Run script: %(script)s" % self.run
        self.headerText  = "Xmipp Protocol: %s\n" % protDict[self.run['protocol_name']].title
        self.headerText += "Project: %s" % os.path.basename(self.project.projectDir) 
        
    def fillHeader(self):
        self.setTitleAndHeader()
        self.master.title(self.titleText)
        self.fonts = {}
        self.fonts['header'] = tkFont.Font(family=FontName, size=FontSize+2, weight=tkFont.BOLD)
        self.l1 = tk.Label(self.frame, text=self.headerText, fg=SectionTextColor, bg=BgColor, 
                        font=self.fonts['header'], pady=5)
        self.l1.configure(wraplength=WrapLenght)
        self.l1.grid(row=self.getRow(), column=0, columnspan=6, sticky='ew')
        self.citerow = self.getRow()        
            
    def scroll(self, event):
        if event.num == 5 or event.delta < 0:
            count = 1
        if event.num == 4 or event.delta > 0:
            count = -1
        self.canvas.yview("scroll", count, "units")
        
    def createWidgets(self):
        ''' Create the widgets after the protocol header script have been 
        parsed and the variables were created '''
        self.checkSpecialCases()
        citeslist = []
        maintext = ''
        
        for s in self.parser.sections:
            section = self.createSectionWidget(s)
            for var in s.childs:
                if var.hasTag('cite'):
                    citeslist.append(var.getValue())
                self.createWidget(section, var)
        
        #Show usage if found
        if self.parser.hasVariable('Usage'):
            maintext += self.getVarValue('Usage')
        #Show citations if found
        if len(citeslist) > 0:
            maintext += "If you publish results obtained with this protocol, please cite:\n"
            maintext += ' '.join(citeslist).strip()
            
        if len(maintext):
            self.fonts['cites'] = tkFont.Font(family=FontName, size=FontSize-2, weight=tkFont.BOLD)
            label = tk.Label(self.frame, text=maintext, fg=CitationTextColor, bg=BgColor, justify=tk.LEFT,
                          font=self.fonts['cites'], wraplength=WrapLenght)
            label.grid(row=self.citerow, column=0, columnspan=5, sticky='ew')
            
    #-------------------------------------------------------------------
    # GUI Events handling
    #-------------------------------------------------------------------           
        
    def close(self, e=None):
        self.master.destroy()
    
    def toggleExpertMode(self, event=""):
        self.expert_mode = not self.expert_mode
        
        if self.expert_mode:
            text = "Less options"
            strValue = 'True'
        else:
            text = "More options"
            strValue = 'False'
        self.setVarValue('ShowExpertOptions', strValue)
        self.btnToggleExpert.config(text=text)
        self.checkVisibility()
        
    def save(self, event=""):
        from protlib_sql import IntegrityError
         
        if not self.validateInput():
            return False
        try:
            runName = self.getVarValue('RunName')
            # Validate runName only contains letters, numbers or underscore:
            import re
            if re.match('\w+$', runName) is None:
                self.showError("Invalid Run Name", 'RUN NAME should only contains letters, numbers or underscore')
                return False
            
            if runName != self.inRunName:
                self.run['run_name'] = runName
                self.run['script'] = self.project.getRunScriptFileName(self.run['protocol_name'], runName)
            
            
            # update database
            if self.run['script'] != self.run['source']:
                self.project.projectDb.insertRun(self.run)
                self.run['source'] = self.run['script']
                self.inRunName = runName
            else:
                if self.run['run_state'] in [SqliteDb.RUN_STARTED, SqliteDb.RUN_LAUNCHED] and not self.visualize_mode:
                    self.showError('Save not allowed', 
                              "This run appears to be <RUNNING> or <LAUNCHED>, so you can't save it")
                    return False 
                if not self.visualize_mode:
                    self.run['run_state'] = SqliteDb.RUN_SAVED
                self.project.projectDb.updateRun(self.run)
            # write header file if not error
            self.parser.save(self.run['script'])
    
            if self.saveCallback:
                self.saveCallback()
                
            FlashMessage(self.master, 'Saved successfully.', delay=1)
        except IntegrityError, ie: 
            if 'columns run_name, protocol_name are not unique' in str(ie):
                self.showError("Error saving run parameters", 
                               "The RUN name <%s> \nalready exists in the PROJECT." % getExtendedRunName(self.run))
            else:
                self.showError('Database error', str(ie))
        except Exception, e:
            self.showError("Error saving run parameters", str(e))
            raise e
        return True
    
    def showError(self, title, e):
        showError(title, str(e), parent=self.master)
    
    def validateInput(self):
        errors = []
        for s in self.sectionslist:
            expanded = s.isSectionExpanded()
            if expanded:   
                for w in s.childwidgets:
                    errors += w.validate()
        if len(errors) > 0:
            showError("Validation ERRORS", '\n'.join(errors), parent=self.master)
            return False
        return True
        
    def validateProtocol(self, prot):
        self.validate_errors = []
        def validateBase():
            self.validate_errors = prot.validateBase()
        
        FlashMessage(self.master, 'Validating data...', func=validateBase)
                
        if len(self.validate_errors) > 0:
            showError("Validation ERRORS", '\n'.join(self.validate_errors), parent=self.master)
            return False
        return True
    
    def getProtocol(self):
        prot = self.project.getProtocolFromModule(self.run['script'])
        prot.parser = self.parser
        return prot
    
    def saveExecute(self, event=""):
        if not self.save(): #Validate user input
            return
        prot = self.getProtocol()        
        if self.visualize_mode:
            prot.master = self.master
            prot.visualize()
        else:
            if self.validateProtocol(prot):
                doRun = True
                warnings = prot.warningsBase()
                if len(warnings) > 0:
                    msg = "There are some warnings:\n%s\n\nDo you want to proceed?" % '\n'.join(warnings)
                    doRun = askYesNo("Confirm execution", msg, self.master)
                if doRun:
                    args = 'xmipp_python %s --no_confirm &' % self.run['script']
                    Popen(args, shell=True)
                    self.master.destroy()
    
    def selectFromList(self, var, list, callback=None):
        from protlib_wizard import wizardSelectFromList
        selected=wizardSelectFromList(self.master,self.frame, list)
        if selected is not None:
            var.setTkValue(selected)
            if callback is not None:
                callback(self)
        
    def showHelp(self, helpmsg):
        showInfo("Help", helpmsg, parent=self.master)
    
    def getRow(self):
        row = self.lastrow
        self.lastrow += 1
        return row
    
    def changeFont(self, event=""):
        changeFontSize(self.style.Font, event, self.style.MinFontSize, self.style.MaxFontSize)
        centerWindows(self.master, self.resize() )
     
    def checkProgramAvailable(self, varName, program):
        if self.hasVar(varName):    
            var = self.parser.getVariable(varName) 
            if which(program) == '':
                var.setValue('False')
                var.setTag('hidden')
                
    def checkSpecialCases(self):
        ''' check special variables which values and visibility can 
        change depending on the conditions '''
        if self.hasVar('Behavior'):
            var = self.parser.getVariable('Behavior')
            var.setValue('Resume')
            if not exists(self.run['script']) or not exists(self.getProtocol().WorkingDir):
                var.setTag('hidden')
        
        from protlib_utils import loadLaunchModule
        launch = loadLaunchModule()
        # Check if mpi and queue are availables
        self.checkProgramAvailable('NumberOfMpi', launch.MpiProgram)
        self.checkProgramAvailable('SubmitToQueue', launch.Program)
        
    def checkVisibility(self, *args):
        for s in self.sectionslist:
            expanded = not s.var.hasQuestion() or self.expandCollapseSection(s)
            # Check the section visibility and check childs if needed                
            if s.checkVisibility() and expanded:
                visible_child = False #Check there is at least one visible child 
                for w in s.childwidgets:
                    if w.checkVisibility():
                        visible_child = True
                if not visible_child:
                    s.display(False)
        #self.resize()
        self.updateScrollRegion() 

    def fillButtons(self):
        row = self.getRow()
        row += 3
        self.btnToggleExpert = self.addButton("More options", self.toggleExpertMode, 5, row, 1, 'ew')
        self.addButton("Save", self.save, 0, row, 2, 'w')
        self.btnExecute = self.addButton("Save & Execute", self.saveExecute, 7, row, 3, 'w')
        self.addButton("Close", self.close, 0, row, 4, 'w')
        
    def addBindings(self):
        self.master.bind('<Alt_L><c>', self.close)
        self.master.bind('<Alt_L><o>', self.toggleExpertMode)
        self.master.bind('<Alt_L><s>', self.save)
        self.master.bind('<Alt_L><e>', self.saveExecute)
        self.master.bind('<Alt_L><r>', self.saveExecute)
        self.master.bind('<Alt_L><plus>', self.changeFont)
        self.master.bind('<Alt_L><minus>', self.changeFont)
        self.master.bind('<Return>', self.checkVisibility)
        # with Windows OS
        self.master.bind("<MouseWheel>", self.scroll)
        # with Linux OS
        self.master.bind("<Button-4>", self.scroll)
        self.master.bind("<Button-5>", self.scroll)
        
    def hasVar(self, varName):
        return self.parser.hasVariable(varName)
    
    def getVarValue(self, varName):
        ''' Return the value of a variable give the name'''
        if self.hasVar(varName):
            var = self.parser.getVariable(varName)
            var.updateValue()
            return var.getValue()
            #return self.parser.getValue(varName)
        return None
    
    def getVarlistValue(self, varList):
        ''' Return a list of values of given variables names'''
        return [self.getVarValue(v) for v in varList]
    
    def getVarLiteralValue(self, varName):
        if self.hasVar(varName):
            var = self.parser.getVariable(varName)
            return var.getLiteralValue()
        return None
    
    def setVarValue(self, varName, value):
        ''' Set the value of some var'''
        if self.hasVar(varName):
            self.parser.getVariable(varName).setTkValue(value)
            
    def setVarlistValue(self, varList, valueList):
        ''' Return a list of values of given variables names'''
        for var, value in zip(varList, valueList):
            self.setVarValue(var, value)
        
    def createGUI(self, project, run, master=None, saveCallback=None, visualize_mode=False):
        try:
            self.createBasicGUI(master)
            self.initialize(project, run, saveCallback)
            self.createScrollableCanvas()
            self.fillHeader()
            self.parser = ProtocolParser(self.run['source'])
            self.createWidgets()
            self.master.update_idletasks()
            maxWidth = max([s.frame.winfo_width() for s in self.sectionslist])
            #self.maxLabelWidth = 300
            for section in self.sectionslist:
                section.frame.grid_columnconfigure(0, minsize=maxWidth)
                section.content.grid_columnconfigure(0, minsize=self.maxLabelWidth)
                #section.content.grid_columnconfigure(1)
                #section.content.grid_columnconfigure(2, weight=1)
            #Set the run_name
            if self.hasVar('RunName'):
                self.inRunName = run['run_name']
                self.setVarValue('RunName', self.inRunName)
                
            self.visualize_mode = visualize_mode
            self.fillGUI()
        except Exception, e:
            errMsg = "Couldn't create GUI. ERROR: %s\n" % e
            showError("GUI Creation Error", errMsg,self.master)
            self.Error = e
            raise
        
    def fillGUI(self):        
        if self.visualize_mode and not self.hasVisualizeOptions:
            return
        self.fillButtons()
        self.addBindings()    
        self.master.update_idletasks()  
        self.expert_mode = self.hasVar('ShowExpertOptions') and self.getVarValue('ShowExpertOptions') == 'True'
        self.checkVisibility()
        
    def launchGUI(self):
        if self.Error:
            return
        if self.visualize_mode and not self.hasVisualizeOptions:
            self.getProtocol().visualize()
        else:
            BasicGUI.launchGUI(self)
    
# Class based on ProtocolGUI but special for single 
# programs execution
class ProgramGUI(ProtocolGUI):
    def __init__(self, script):
        self.script = script
        self.loadModule()
        self.msgWin = None
    
    def loadModule(self):
        scriptMod = loadModule(self.script)
        for k in dir(scriptMod):
            if not (k.startswith('_') and not k.startswith('_K_')):
                obj = getattr(scriptMod, k)
                if not callable(obj):
                    setattr(self, k, obj)
                
    def workingDirPath(self, value):
        return value
    
    def save(self):
        self.parser.save(self.script)
        self.loadModule()

    def getCmd(self):
        from protocol_program import buildCommandLine
        self.cmdStr = "%s %s" % (self.ProgramName, buildCommandLine(self))
        
    def closeMsgWin(self, e=None):
        self.msgWin.destroy()
        self.msgWin = None
        
    def showMsg(self, title, msg):
        '''Open message window if not exist to show a message '''
        if not self.msgWin: # Create the msg windows if not exists
            self.msgWin = tk.Toplevel(self.master)
            self.msgText = TaggedText(self.msgWin, width=70, height=7, background="white")
            self.msgText.frame.grid(row=0, column=0, columnspan=3, padx=5)
            XmippButton(self.msgWin, text="Close", width=6, command=self.closeMsgWin).grid(row=1, column=2, pady=5, padx=5, sticky='e')
            centerWindows(self.msgWin, refWindows=self.master)
        
        self.msgWin.title(title)
        self.msgText.clear()
        self.msgText.addText(msg)
        self.msgWin.lift() # raise to the top
        
    def showCmd(self, e=None):
        self.save()
        self.getCmd()
        self.showMsg('Command line', self.cmdStr)
        
    def validate(self):
        from subprocess import PIPE
        cmd = self.cmdStr + ' --xmipp_validate_params'
        ps = Popen(cmd, shell=True, stderr=PIPE)
        return ps.communicate()[1] # return output of stderr

    def saveExecute(self, e=None):
        self.save()       
        self.getCmd()
        errors = self.validate()
        if len(errors):
            self.showMsg("Errors", errors)
        else:
            self.master.destroy()
            print "Running: ", greenStr(self.cmdStr)
            os.system(self.cmdStr)    
            os.remove(self.script + 'c')#.pyc
            os.remove(self.script) #remove script    
        
    def setTitleAndHeader(self):
        self.titleText = 'Program: %s' % self.ProgramName
        self.headerText = self.titleText
        
    def fillButtons(self):
        row = self.getRow()
        row += 3
        self.btnToggleExpert = self.addButton("More options", self.toggleExpertMode, 5, row, 1, 'ew')
        self.addButton("Show command", self.showCmd, 0, row, 2, 'w')
        self.btnExecute = self.addButton("Execute", self.saveExecute, 0, row, 3, 'w')
        self.addButton("Close", self.close, 0, row, 4, 'w')
        
    def addBindings(self):
        ProtocolGUI.addBindings(self)
        self.master.bind('<Alt_L><s>', self.showCmd)
        
def launchProgramGUI(script):
    run = {'source': script}
    root = tk.Tk()
    gui = ProgramGUI(script)
    gui.createGUI(None, run, root)
    gui.launchGUI()
