#!/usr/bin/env python
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
import Tkinter as tk
import tkMessageBox
import tkFont
from protlib_base import protocolMain, getProtocolFromModule
from protlib_utils import loadModule
from protlib_gui_ext import centerWindows
from protlib_filesystem import getXmippPath
from config_protocols import protDict
from config_protocols import FontName, FontSize, MaxHeight, MaxWidth, WrapLenght
from config_protocols import LabelTextColor, SectionTextColor, CitationTextColor
from config_protocols import BgColor, EntryBgColor, SectionBgColor, LabelBgColor, ButtonActiveBgColor                         

class ProtocolStyle():
    ''' Class to define some style settings like font, colors, etc '''
    def __init__(self, configModuleName=None):        
        #Font
        self.FontName = "Helvetica"
        self.FontSize = 10
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

class AutoScrollbar(tk.Scrollbar):
    '''A scrollbar that hides itself if it's not needed.'''
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            self.tk.call("grid", "remove", self)
        else:
            self.grid()
        tk.Scrollbar.set(self, lo, hi)

def createSection(parent, text):
    frame = tk.Frame(parent, bd=2, relief=tk.RAISED, bg=SectionBgColor)
    frame.columnconfigure(0, weight=1)
    label = tk.Label(frame, text=text, fg=LabelTextColor, bg=SectionBgColor)
    label.grid(row=0, column=0, sticky=tk.W)        
    content = tk.Frame(frame, bg=LabelBgColor, bd=0)
    content.grid(row=1, column=0, columnspan=5, sticky=tk.NSEW, ipadx=5, ipady=5)  
    return (frame, label, content)     

class BasicGUI(): 
    def __init__(self):
        pass
        
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
        
    def resize(self):
        height = min(self.frame.winfo_reqheight() + 25, MaxHeight)
        width = min(self.frame.winfo_reqwidth() + 25, MaxWidth)
        x = self.frame.winfo_x()
        y = self.frame.winfo_y()
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
                      bg=self.style.BgColor)
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
                          selectcolor=self.style.ButtonActiveBgColor)
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
        
sepLine = "#------------------------------------------------------------------------------------------\n" 
        
class ProtocolVariable():
    '''Store information about Protocols variables.'''
    def __init__(self):
        self.name = None
        self.value = None
        self.tkvar = None
        self.comment = None
        self.help = None
        self.tags = {}
        self.conditions = {}
        self.is_string = False    
        self.tktext = None #special text widget                       
    
    def setTags(self, tags):
        for k, v in tags:
            self.tags[k] = v
            
    def isExpert(self):
        return 'expert' in self.tags.keys()
    
    def isVisualize(self):
        return 'visualize' in self.tags.keys()
    
    def isSection(self):
        return 'section' in self.tags.keys()
    
    def isHidden(self):
        return 'hidden' in self.tags.keys()
    
    def getValue(self):
        if self.tktext:
            return self.tktext.get(1.0, tk.END)
        return self.tkvar.get() 

    def setValue(self, value):
        self.tkvar.set(value)
    
    def getLines(self):   
        if self.isSection():
            lines = [sepLine + self.commentline + '\n' + sepLine]
        else:    
            lines = [self.commentline + '\n']
        if 'text' in self.tags.keys(): #Store the value in the comment field
            self.help = self.getValue()
        if self.help:
            lines.append('"""%s"""\n' %  self.help)
        if self.is_string:
            template = '%s = "%s"\n\n'
        else:
            template = '%s = %s\n\n'  
        if self.tkvar:   
            lines.append(template % (self.name, self.getValue()))
        return lines
    
class ProtocolWidget():
    def __init__(self, master, var):
        self.master = master
        self.widgetslist = [] # Store real tk widgets
        self.childwidgets = [] # This will only be used for sections
        self.variable = var
       
    def satisfiesCondition(self):
        if not (self.variable and self.variable.conditions):
            return True
        for k, v in self.variable.conditions.iteritems():
            var = self.master.variablesDict[k]
            if var.getValue() != v:
                return False
        return True
     
    def checkVisibility(self):
        show = self.variable.isVisualize() == self.master.visualize_mode
        if show and self.variable.isExpert():
            show = self.master.expert_mode
        if show:
            show = self.satisfiesCondition() #has condition, check it
        self.display(show)    

            
    def display(self, value):
        for w in self.widgetslist:
            if value:
                w.grid()
            else:
                w.grid_remove()

""" Guidelines for python script header formatting:

The idea of the GUI class is that it provides a graphical editor of
only the variables contained in the header of a python script, and
that it can also be used to launch the script.

Usage from main:
    gui = ProtocolGUI()
    gui.createGUI(script)
    gui.fillGUI()
    gui.launchGUI()
    
where script contains the filename of the python script that will be
parsed. Or if you already has a GUI environmet, you can launch the
protocol GUI in another  window by simply creating a TopLevel windows
and passing as root to createGUI function. Following an example
    top = TopLevel()
    gui = ProtocolGUI()
    gui.createGUI(script, top)
    gui.fillGUI()
    gui.launchGUI()
    
Where the header of script.py should be organized as follows:

Obligatory:
    * Include a {begin_of_header} label at the begining of the header
      and {end_of_header} label at the end of the header
    * Variable declaration (visible in GUI): Variablename = XXX
          o If XXX is True or False, the variable is considered as a
            Boolean (yes/no button in GUI)
          o If XXX starts with " or ', the variable is considered as
            a string (text field in GUI) or a list option (see {list} tag}
          o If XXX is a number, the variable is considered as a number
            (text field in GUI) 
    * The first single comment line above each variable (starting with a #)
      will be displayed in the GUI, as a label left from the variable entry field
    * More detailed help for each variable can be placed between the comment
      line and the variable declaration line using \""" ...\""".
      This text will become visible in the GUI by pressing a HELP
      button, at the right side of the variable entry field.
Optional:
    * There are extra tags in the comment that have some semantic for the GUI:
        - {section} defines a section with the comment, not a variable      
        - {has_question} this tag can be used with a section and take the next 
              boolean variables as a checkbutton for show/hide the entire section box
        - {expert} mark this variable or section as "expert", not shown in the GUI
              by default, showed when pressed the "Show expert options" button
        - {file} or {dir} marks the option as a Filename or Directory, and will add 
              a corresponding "Browse" button to that option
        - {view} this will add an extra button to view/explore some files or folders
        - {run}(prot1, prot2) Select as input previous runs from other(or same) protocol
        - {text} Will display a text box with more space for writing
        - {list}(option A, option B, option C) marks the option as a radio-list button. 
              The selected variable should be one of the options indicated.
        - {condition}(option=True)
          {condition}(selection=optionA) marks an option or section dependent of some
              condition, it can be used with boolean options or for list-selections
        - {please_cite} will display
              a message at the top of the protocol stating -If you publish results obtained with
              this protocol, please cite-, followed by the text on rest of the comment line.
              If more than one citation lines are present, they will all be displayed.
              DONT use very long citations, as this will results in an ugly gui.
        - {hidden} label on the comment line (#) marks the option as -hidden-
        - {wizzard} (WizzardFunction) this will serve to plugin a graphical wizzard
             to select some parameter
"""       

class ProtocolGUI(BasicGUI):
    def init(self):
        self.variablesDict = {}       
        self.pre_header_lines = []
        self.header_lines = []
        self.post_header_lines = []
        self.expert_mode = True
        self.visualize_mode = False
        self.lastrow = 0
        self.widgetslist = []
        self.sectionslist = [] # List of all sections
        self.citeslist = []
        # Script title
        self.programname = os.path.basename(self.run['source'].replace('.py', ''))
        self.maxLabelWidth = 0

    #-------------------------------------------------------------------
    # Widgets creation and GUI building
    #-------------------------------------------------------------------  
        
    def addButton(self, text, cmd, underline, row, col, sticky, imageFilename=None, parent=None, tip=None):
        f = tkFont.Font(family=self.style.FontName, size=self.style.ButtonFontSize, weight=tkFont.BOLD)
        helpImage = None
        if parent is None:
            parent = self.frame
        if imageFilename:
            try:
                imgPath = os.path.join(getXmippPath('resources'), imageFilename)
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
        rb = tk.Radiobutton(parent, text=text, variable=var.tkvar, value=value, bg=self.style.LabelBgColor)#, command=self.checkVisibility)
        rb.grid(row=row, column=col, sticky='w')
        w.widgetslist.append(rb)
        return rb
        
    def createSectionWidget(self, w, var):
        w.frame, w.label, w.content = createSection(self.frame, var.comment)
        w.name = var.comment
        w.widgetslist.append(w.frame)
        
    def expandCollapseSection(self, section):
        if section.tkvar.get() == 'False':
            section.content.grid_remove()
        else:
            section.content.grid(row=1, column=0, columnspan=5, sticky='nsew')
        self.updateScrollRegion()
            
    def createWidget(self, var):
        w = ProtocolWidget(self, var)  
        self.widgetslist.append(w)  
        label_row = row = self.getRow() # Get row where to place the widget        
        label_text = var.comment
        label_color = self.style.LabelTextColor
        label_bgcolor = self.style.LabelBgColor        
        
        if 'section' in var.tags.keys():
            self.createSectionWidget(w, var)
            w.frame.grid(row=label_row, column=0, columnspan=5, 
                         sticky='ew', pady=5, padx=(10,5))
            w.has_question = 'has_question' in var.tags.keys()
            self.lastSection = w
            section = w
            self.sectionslist.append(w)
        else:
            section = self.lastSection
            section.childwidgets.append(w)
            frame = section.content
            #widgets inherit expert from section and its conditions 
            if section.variable.isExpert():
                var.tags['expert'] = section.variable.tags['expert']
            if section.variable.isVisualize():
                var.tags['visualize'] = section.variable.tags['visualize']
            for k, v in section.variable.conditions.iteritems():
                var.conditions[k] = v
                
        keys = var.tags.keys()
                
        if 'expert' in keys:
            label_bgcolor = self.style.ExpertLabelBgColor
           
        if 'condition' in keys:
            conditions = var.tags['condition'].split(',')
            for cond in conditions:
                cond_name, cond_value = cond.split('=')
                var.conditions[cond_name.strip()] = cond_value.strip()
         
        if 'text' in keys:
            #scrollbar = AutoScrollbar(frame)
            #scrollbar.grid(row=label_row, column=6, sticky=NS)
            var.tktext = tk.Text(frame, width=66, height=10, wrap=tk.WORD, bg=EntryBgColor)#, yscrollcommand=scrollbar.set, bg=EntryBgColor)
            #scrollbar.config(command=var.tktext.yview)
            var.tktext.grid(row=label_row, column=0, columnspan=5, sticky='ew', padx=(10, 0), pady=(10, 0))
            var.tktext.insert(tk.END, var.help)
            w.widgetslist.append(var.tktext)
            #w.widgetslist.append(scrollbar)
            return w
        # Treat variables
        if var.value:
        #Escape string literals
            var_column = 1
            var.tkvar = tk.StringVar()   
            
            if var.value.startswith('"') or var.value.startswith("'"):
                var.value = var.value.replace('"', '')
                var.value = var.value.replace("'", '')
                var.is_string = True
                
            var.tkvar.set(var.value)
            
            if 'hidden' in keys:
                return w
            
            if var.value == 'True' or var.value == 'False':
                #Check if that boolean variable is the question of the section
                if section.has_question and len(section.childwidgets) == 1:
                    section.tkvar = var.tkvar
                    #Label(section.frame, text=label_text, bg=SectionBgColor).grid(row=0, column=1, padx=(5, 0))
                    chb = tk.Checkbutton(section.frame, text=label_text, variable=var.tkvar, 
                                onvalue='True', offvalue='False',
                                command=lambda:self.expandCollapseSection(section),
                                bg=SectionBgColor, activebackground=ButtonActiveBgColor)
                    chb.grid(row=0, column=1, padx=(5, 0))
                    self.expandCollapseSection(section)
                    return w
                else:
                    self.addRadioButton(w, var, 'Yes', 'True', row, var_column, frame)  
                    self.addRadioButton(w, var, 'No', 'False', row, var_column + 1, frame) 
            elif 'list' in keys:
                opts = var.tags['list'].split(',')
                for o in opts:
                    o = o.strip()
                    self.addRadioButton(w, var, o, o, row, var_column, frame)
                    row = self.getRow()
            elif 'text' in keys:
                scrollbar = tk.Scrollbar(frame)
                scrollbar.grid(row=label_row+1, column=1, sticky='ns')
                var.tktext = tk.Text(frame, width=66, height=10, wrap=tk.WORD, yscrollcommand=scrollbar.set, bg=EntryBgColor)
                scrollbar.config(command=var.tktext.yview)
                var.tktext.grid(row=label_row+1, column=0, columnspan=5, sticky='ew', padx=(10, 0))
                w.widgetslist.append(var.tktext)
                var.tktext.insert(tk.END, var.value)            
            else: #Add a text Entry
                entry = tk.Entry(frame, textvariable=var.tkvar, bg=EntryBgColor)
                entry.grid(row=row, column=var_column, columnspan=2, sticky='ew')
                w.widgetslist.append(entry)
                args = None
                if 'file' in keys:
                    args = ['Browse', lambda: self.browse(var.tkvar, True), 'fileopen.gif', 'Browse file']
                elif 'dir' in keys:
                    args = ['Browse', lambda: self.browse(var.tkvar, False), 'folderopen.gif', 'Browse folder']
                elif 'run' in keys:
                    protocols = var.tags['run'].split(',')
                    runs = []
                    for p in protocols:
                        runs += self.project.projectDb.selectRunsByProtocol(p)
                    list = ["%s_%s" % (r[5], r[1]) for r in runs]
                    args = ['Select Run', lambda: self.selectFromList(var, list), 'wizzard.gif', 'Select run']
                elif 'blocks' in keys:
                    #md = self.variablesDict[var.tags['blocks']].getValue()
                    args = ['Select Blocks', lambda: self.selectFromList(var, ['block1', 'block2', 'block3']), 'wizzard.gif', 'Select blocks']
                
                if args:
                    btn = self.addButton(args[0], args[1], -1, label_row, var_column+2, 'nw', args[2], frame, args[3])
                    w.widgetslist.append(btn)
                
                if 'view' in keys:
                    btn = self.addButton("View", lambda: self.viewFiles(), -1, label_row, var_column+3, 'nw', 'visualize.gif', frame, 'View file')
            if var.help:
                btn = self.addButton("Help", lambda: self.showHelp(var.help.replace('"', '')), -1, label_row, var_column+4, 'nw', 'help.gif', frame, 'Show info')
                w.widgetslist.append(btn)
            if var.name == 'RunName':
                label_text += ' %s_' % self.run['protocol_name']
                    
            label = tk.Label(frame, text=label_text, fg=label_color, bg=label_bgcolor)
            label.grid(row=label_row, column=0, sticky='e', padx=(5, 10))
            self.maxLabelWidth = max(self.maxLabelWidth, label.winfo_reqwidth())
            w.widgetslist.append(label)
        
        return w
        
    def fillHeader(self):
        self.master.title("Run script: %(script)s" % self.run)
        headertext  = "Xmipp Protocol: %s\n" % protDict.protocolDict[self.run['protocol_name']].title
        headertext += "Project: %s" % self.project.projectDir 
        self.fonts = {}
        self.fonts['header'] = tkFont.Font(family=FontName, size=FontSize+2, weight=tkFont.BOLD)
        self.l1 = tk.Label(self.frame, text=headertext, fg=SectionTextColor, bg=BgColor, 
                        font=self.fonts['header'], pady=5)
        self.l1.configure(wraplength=WrapLenght)
        self.l1.grid(row=self.getRow(), column=0, columnspan=6, sticky='ew')
        self.citerow = self.getRow()        
            
    def browse(self, var, isFile):
        import tkFileDialog
        if isFile:
            filename = tkFileDialog.askopenfilename(title="Choose file", parent=self.master)
        else:
            filename = tkFileDialog.askdirectory(title="Choose directory", parent=self.master)
        if len(filename) > 0:
            var.set(os.path.relpath(filename))
         
    def scroll(self, event):
        if event.num == 5 or event.delta == -120:
            count = 1
        if event.num == 4 or event.delta == 120:
            count = -1
        self.canvas.yview("scroll", count, "units")
    #-------------------------------------------------------------------
    # Reading and parsing script
    #-------------------------------------------------------------------
    def readProtocolScript(self):
        begin_of_header = False
        end_of_header = False    
        script = self.run['source']    
        f = open(script, 'r')
        for line in f:
            #print "LINE: ", line
            if not begin_of_header:
                self.pre_header_lines.append(line)
            elif not end_of_header:
                #print "LINE: ", line
                self.header_lines.append(line)
            else:
                self.post_header_lines.append(line)                
            if line.find('{begin_of_header}') != -1:
                begin_of_header = True
            if line.find('{end_of_header}') != -1:
                end_of_header = True
        f.close()
        
        if not begin_of_header:
            raise Exception('{begin_of_header} tag not found in protocol script: %s' % script)
        if not end_of_header:
            raise Exception('{end_of_header} tag not found in protocol script: %s' % script)
                
    def parseHeader(self):
        #REGURLAR EXPRESSION TO PARSE VARIABLES DEFINITION
        import re
        #Comment regex, match lines starting by # and followed by tags with values
        #of the form {tag}(value) and ending with the comment for the GUI    
        reComment = re.compile('#\s*((?:{\s*\w+\s*}\s*(?:\([^)]*\))?)*)?\s*(.*)')
        #This take all tags and values from previous one
        reTags = re.compile('(?:{\s*(\w+)\s*}\s*(?:\(([^)]*)\))?\s*)')
        #This is the regular expression of a Variable
        #possible values are: True, False, String with single and double quotes and a number(int or float) 
        #reVariable = re.compile('(\w+)\s*=\s*(True|False|".*"|\'.*\'|\d+|)')
        reVariable = re.compile('(\w+)\s*=\s*(True|False|".*"|\'.*\'|[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?|)')
        
        self.variablesDict = {}
        self.widgetslist = []
        index = 0;
        count = len(self.header_lines)
        while index < count:
            line = self.header_lines[index].strip()
            index += 1
            match = reComment.search(line) #Parse the comment line
            if match:
                v = ProtocolVariable()
                v.comment = match.group(2)
                v.commentline = line
                
                if match.group(1) != '':
                    tags = reTags.findall(match.group(1))
                    v.setTags(tags)
                    
                is_section = v.isSection()
                if not is_section:
                    #This is a variable, try to get help string
                    helpStr = ''
                    countQuotes = 0 # count starting and endi
                    if index < count and self.header_lines[index].strip().startswith('"""'):
                        while index < count and countQuotes < 2:
                            line = self.header_lines[index].strip()
                            if line.startswith('"""'):
                                countQuotes += 1
                            if len(line) > 4 and line.endswith('"""'):
                                countQuotes += 1
                            line = line.replace('"""', '')
                            if len(line) > 0:
                                helpStr += line + '\n'
                            index += 1
                        v.help = helpStr
                        
                    if 'please_cite' in v.tags.keys():
                        self.citeslist.append(v.help)
                    else:
                        if index < count:
                            line = self.header_lines[index].strip()
                            match2 = reVariable.match(line)
                            if match2:
                                v.name, v.value = (match2.group(1).strip(), match2.group(2).strip())
                                #print "DEBUG_JM: v.name: '%s', v.value: '%s'" % (v.name, v.value)
                                self.variablesDict[v.name] = v
                                index += 1
                if is_section or v.name or 'text' in v.tags.keys():
                    self.createWidget(v)
        #Update if citations found
        if len(self.citeslist) > 0:
            citetext = "If you publish results obtained with this protocol, please cite:\n"
            citetext += '\n'.join(self.citeslist)
            self.fonts['cites'] = tkFont.Font(family=FontName, size=FontSize-2, weight=tkFont.BOLD)
            label = tk.Label(self.frame, text=citetext, fg=CitationTextColor, bg=BgColor,
                          font=self.fonts['cites'], wraplength=WrapLenght)
            label.grid(row=self.citerow, column=0, columnspan=5, sticky='ew')
            
    #-------------------------------------------------------------------
    # GUI Events handling
    #-------------------------------------------------------------------           
        
    def close(self, event=""):
        self.master.destroy()
    
    def toggleExpertMode(self, event=""):
        self.expert_mode = not self.expert_mode
        
        if self.expert_mode:
            text = "Hide Expert Options"
            strValue = 'True'
        else:
            text = "Show Expert Options"
            strValue = 'False'
        if 'ShowExpertOptions' in self.variablesDict.keys():
            self.variablesDict['ShowExpertOptions'].setValue(strValue)
        self.btnToggleExpert.config(text=text)
        self.checkVisibility()
        
    def toggleVisualizeMode(self, event=''):
        self.visualize_mode = not self.visualize_mode
        if self.visualize_mode:
            text = "Show Run Options"
        else:
            text = "Show Visualize Options"
        self.btnToggleVisualize.config(text=text)    
        self.checkVisibility()
    
    def save(self, event=""):
        try:
            runName = self.getRunName()
            if runName != self.inRunName:
                self.run['run_name'] = runName
                self.run['script'] = self.project.getRunScriptFileName(self.run['protocol_name'], runName)
            #update database
            if self.run['script'] != self.run['source']:
                self.project.projectDb.insertRun(self.run)
                self.run['source'] = self.run['script']
                self.inRunName = runName
            else:
                self.project.projectDb.updateRun(self.run)
    
            #print "* Saving script: %s" % self.run['script']
            f = open(self.run['script'], 'w')
            #f = sys.stdout
            f.writelines(self.pre_header_lines)
            if len(self.citeslist) > 0:
                f.writelines('#{please_cite}\n"""\n')
                f.writelines(self.citeslist)
                f.writelines('"""\n')
            for w in self.widgetslist:
                wlines = w.variable.getLines()
                f.writelines(wlines)
            f.writelines(sepLine + '# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE\n') 
            f.writelines(self.post_header_lines)
            f.close()
            os.chmod(self.run['script'], 0755)
        except Exception, e:
            tkMessageBox.showerror("Error saving run parameters", str(e), parent=self.master)
        if self.saveCallback:
            self.saveCallback()
    
    def confirmDeleteWorkingDir(self):
        if 'DoDeleteWorkingDir' in self.variablesDict and self.variablesDict['DoDeleteWorkingDir'].getValue() == 'True':
            return tkMessageBox.askyesno("Confirm DELETE", "Working dir '%s' will be DELETED. Do you want to continue?" % self.variablesDict['WorkingDir'].getValue())
        return True
    
    def validateProtocol(self):
        prot = getProtocolFromModule(self.run['script'], self.project)
        errors = prot.validateBase()        
        if len(errors) > 0:
            tkMessageBox.showerror("Validation ERRORS", '\n'.join(errors), parent=self.master)
            return None
        return prot
    
    def saveExecute(self, event=""):
        self.save() 
        prot = self.validateProtocol()
        if not prot is None:
            warnings=prot.warningsBase()
            if len(warnings)==0 or tkMessageBox.askyesno("Confirm execution",'\n'.join(warnings), parent=self.master):
                os.system('python %s --no_confirm >>%s 2>>%s &' % (self.run['script'], prot.Out, prot.Err) )
                self.master.destroy() 
    
    def viewFiles(self):
        tkMessageBox.showinfo("Visualize", "This should open ImageJ plugin to display files", parent=self.master)
        
    def selectFromList(self, var, list):
        from protlib_gui_ext import ListboxDialog
        d = ListboxDialog(self.frame, list, selectmode=tk.SINGLE)
        if len(d.result) > 0:
            index = d.result[0]
            var.setValue(list[index])
        
    def showHelp(self, helpmsg):
        tkMessageBox.showinfo("Help", helpmsg, parent=self.master)
    
    def getRow(self):
        row = self.lastrow
        self.lastrow += 1
        return row
    
    def changeFont(self, event=""):
        deltha = 2
        if event.char == '-':
            deltha = -2
        size = self.style.Font['size']
        new_size = size + deltha
        if new_size >= self.style.MinFontSize and new_size <= self.style.MaxFontSize:
            self.style.Font.configure(size=new_size)
        centerWindows(self.master, self.resize() )
        
    def checkVisibility(self, event=""):
        for s in self.sectionslist:
            if s.has_question:
                self.expandCollapseSection(s)
            s.checkVisibility()
            for w in s.childwidgets:
                w.checkVisibility()
        centerWindows(self.master, self.resize() )
        self.updateScrollRegion() 

    
    def fillButtons(self):
        row = self.getRow()
        row += 3
        self.btnToggleVisualize = self.addButton("Show Visualize Options", self.toggleVisualizeMode, 0, row, 0, 'w')
        self.btnToggleExpert = self.addButton("Show Expert Options", self.toggleExpertMode, 12, row, 1, 'ew')
        self.addButton("Save", self.save, 0, row, 3, 'w')
        self.addButton("Save & Execute", self.saveExecute, 7, row, 4, 'w')
        
    def addBindings(self):
        self.master.bind('<Alt_L><c>', self.close)
        self.master.bind('<Alt_L><o>', self.toggleExpertMode)
        self.master.bind('<Alt_L><v>', self.toggleVisualizeMode)
        self.master.bind('<Alt_L><s>', self.save)
        self.master.bind('<Alt_L><e>', self.saveExecute)
        self.master.bind('<Alt_L><r>', self.saveExecute)
        self.master.bind('<Alt_L><plus>', self.changeFont)
        self.master.bind('<Alt_L><minus>', self.changeFont)
        # with Windows OS
        self.master.bind("<MouseWheel>", self.scroll)
        # with Linux OS
        self.master.bind("<Button-4>", self.scroll)
        self.master.bind("<Button-5>", self.scroll)
        
    def getRunName(self):
        return self.variablesDict['RunName'].getValue()
    
    def setRunName(self, value):
        self.variablesDict['RunName'].setValue(value)
        
    def createGUI(self, project, run, master=None, saveCallback=None, visualize_mode=False):
        self.run = run
        self.saveCallback = saveCallback
        self.project = project
        self.init()        
        self.createBasicGUI(master)
        self.master.option_add("*Font", self.style.Font)
        #self.columnspantextlabel = 3
        #self.columntextentry = 3
        
        self.readProtocolScript()
        self.createScrollableCanvas()
        self.fillHeader()
        self.parseHeader()
        self.master.update_idletasks()
        maxWidth = max([s.frame.winfo_width() for s in self.sectionslist])
        #self.maxLabelWidth = 300
        for section in self.sectionslist:
            section.frame.grid_columnconfigure(0, minsize=maxWidth)
            section.content.grid_columnconfigure(0, minsize=self.maxLabelWidth)
            #section.content.grid_columnconfigure(1)
            #section.content.grid_columnconfigure(2, weight=1)
        #Set the run_name
        if self.run:
            self.inRunName = run['run_name']
            self.setRunName(self.inRunName)
            
        self.visualize_mode = visualize_mode
        #self.fillWidgets()
                # Add bottom row buttons
    def fillGUI(self):
        self.fillButtons()
        self.addBindings()    
        self.master.update_idletasks()  
        self.expert_mode = 'ShowExpertOptions'in self.variablesDict and \
                            self.variablesDict['ShowExpertOptions'].getValue() == 'True'
        self.checkVisibility()  
    
