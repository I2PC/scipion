#!/usr/bin/env python
import sys
import os
import string
from Tkinter import *
import tkFont
import tkMessageBox
from protlib_utils import *
from protlib_gui import *

""" Guidelines for python script header formatting:

The idea of the GUI class is that it provides a graphical editor of
only the variables contained in the header of a python script, and
that it can also be used to launch the script.

Usage:
  python xmipp_protocol_gui.py script.py

Where the header of script.py should be organized as follows:

Obligatory:
    * Include a {end_of_header} label at the end of the header
    * Variable declaration (visible in GUI): Variablename=XXX
          o If XXX is True or False, the variable is considered as a
            Boolean (yes/no button in GUI)
          o If XXX starts with \" or \', the variable is considered as
            a string (text field in GUI)
          o If XXX is a number, the variable is considered as a number
            (text field in GUI) 
    * The first single comment line above each variable (starting with a #)
      will be displayed in the GUI, as a label left from the variable entry field

    * More detailed help for each variable can be placed between the comment
      line and the variable declaration line using \""" ...\""".
      This text will become visible in the GUI by pressing a -what's this?-
      button, at the right side of the variable entry field.
      !!!NOTE that the first character in newlines within these comments
      should start with spaces or a tab!!!
    * An {expert} label on the comment line (#) marks the option as -expert-
      (by default not shown in the GUI, unless you press the -show expert options-
      button, at the left side of the variable entry field. Then, they will be
      shown in yellow.
    * A {file} or {dir} label on the comment line (#) marks the option as a
      Filename or Directory, and will add a corresponding Browse-button to the GUI
      Note that this button return absolute paths
    * A {list}|option A|option B|option C| label on the comment line (#) marks 
      the option as a radio-list button. The selected variable should be one of the options indicated.
      The number of different options is not limited.
    * A {hidden} label on the comment line (#) marks the option as -hidden-

Optional:

    * A {please cite} label on a comment line (starting with a #) will display
      a message at the top of the protocol stating -If you publish results obtained with
      this protocol, please cite-, followed by the text on rest of the comment line.
      If more than one citation lines are present, they will all be displayed.
      DONT use very long citations, as this will results in an ugly gui.
    * Include a {section} label on a comment line (starting with a #) to separate
      variables by a blue line + the corresponding title in the GUI 

"""

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
    
    def isSection(self):
        return 'section' in self.tags.keys()
    
    def getValue(self):
        if self.tktext:
            return '""%s""' % self.tktext.get(1.0, END)
        return self.tkvar.get() 
    
    def getLines(self):   
        if self.isSection():
            lines = [sepLine + self.commentline + '\n' + sepLine]
        else:    
            lines = [self.commentline + '\n']
        if self.help:
            lines.append(self.help)
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
        show = True
        c = self.variable.comment
        if self.variable.isExpert():
            show = self.master.expert_mode
        show = show and self.satisfiesCondition() #has condition, check it
        self.display(show)    
            
    def display(self, value):
        for w in self.widgetslist:
            if value:
                w.grid()
            else:
                w.grid_remove()
        for child in self.childwidgets:
            child.display(value)
        
class ProtocolGUI(BasicGUI):
    def init(self, script):
        self.variablesDict = {}       
        self.pre_header_lines = []
        self.header_lines = []
        self.post_header_lines = []
        self.expert_mode = False
        self.scriptname = script
        self.lastrow = 0
        self.widgetslist = []
        self.sectionslist = [] # List of all sections
        self.citeslist = []
        # Script title
        self.programname = os.path.basename(self.scriptname.replace('.py', ''))

    #-------------------------------------------------------------------
    # Widgets creation and GUI building
    #-------------------------------------------------------------------  
    def addSeparator(self, row):
        self.l1 = Label(self.frame, text="", bg=self.style.LabelBgColor)
        self.l1.grid(row=row)
        self.l2 = Frame(self.frame, height=2, bd=1, bg=self.style.SectionTextColor, relief=RIDGE)
        self.l2.grid(row=row + 1, column=0, columnspan=self.columnspantextlabel + 3, sticky=EW)
        self.l3 = Label(self.frame, text="", bg=self.style.LabelBgColor)
        self.l3.grid(row=row + 2)
        self.lastrow += 2
        
    def addButton(self, text, cmd, underline, row, col, sticky):
        f = tkFont.Font(family=self.style.FontName, size=self.style.ButtonFontSize, weight=tkFont.BOLD)
        btn = Button(self.frame, text=text, command=cmd, underline=underline, font=f,
                     bg=self.style.ButtonBgColor, activebackground=self.style.ButtonActiveBgColor)
        btn.grid(row=row, column=col, sticky=sticky)
        return btn
    
    def addRadioButton(self, w, var, text, value, row, col):
        rb = Radiobutton(self.frame, text=text, variable=var.tkvar, value=value, bg=self.style.BgColor, command=self.checkVisibility)
        rb.grid(row=row, column=col, sticky=W)
        w.widgetslist.append(rb)
        return rb
        
    def createWidget(self, var):
        w = ProtocolWidget(self, var)  
        self.widgetslist.append(w)  
        label_row = row = self.getRow() # Get row where to place the widget        
        label_text = var.comment
        label_color = self.style.LabelTextColor
        label_bgcolor = self.style.LabelBgColor        

        if 'section' in var.tags.keys():
            label_text += "\n-----------------------------------------------------------"
            label_color = self.style.SectionTextColor 
            self.lastSection = w
            self.sectionslist.append(w)
        else:
            self.lastSection.childwidgets.append(w)
            #widgets inherit expert from section and its conditions 
            if self.lastSection.variable.isExpert():
                var.tags['expert'] = self.lastSection.variable.tags['expert']
            for k, v in self.lastSection.variable.conditions.iteritems():
                var.conditions[k] = v
                
        keys = var.tags.keys()
                
        if 'expert' in keys:
           label_bgcolor = self.style.ExpertLabelBgColor
           
        if 'condition' in keys:
            conditions = var.tags['condition'].split(',')
            for cond in conditions:
                cond_name, cond_value = cond.split('=')
                var.conditions[cond_name.strip()] = cond_value.strip()
         
        # Treat variables
        if var.value:
        #Escape string literals
            var_column = self.columntextentry
            var.tkvar = StringVar()   
            
            if var.value.startswith('"') or var.value.startswith("'"):
                var.value = var.value.replace('"', '')
                var.value = var.value.replace("'", '')
                var.is_string = True
                
                
            var.tkvar.set(var.value)
            
            if var.value == 'True' or var.value == 'False':
                self.addRadioButton(w, var, 'Yes', 'True', row, var_column)  
                self.addRadioButton(w, var, 'No', 'False', row, var_column + 1)  
            elif 'list' in keys:
                opts = var.tags['list'].split(',')
                for o in opts:
                    o = o.strip()
                    self.addRadioButton(w, var, o, o, row, var_column)
                    row = self.getRow()
            elif 'text' in keys:
                var.tktext = Text(self.frame, width=30, height=10,)
                var.tktext.grid(row=row, column=self.columntextentry, columnspan=2, sticky=W+E)
                w.widgetslist.append(var.tktext)
                var.tktext.insert(END, var.value)
            else: #Add a text Entry
                entry = Entry(self.frame, textvariable=var.tkvar, bg=self.style.EntryBgColor)
                entry.grid(row=row, column=self.columntextentry, columnspan=2, sticky=W+E)
                w.widgetslist.append(entry)
                if 'file' in keys or 'dir' in keys:
                    btn = self.addButton("Browse", lambda: self.browse(var.tkvar, ('file' in keys)), -1, label_row, var_column + 3, NW)
                    w.widgetslist.append(btn)
            if var.help:
                btn = self.addButton("Help", lambda: self.showHelp(var.help), -1, label_row, var_column + 4, NW)
                w.widgetslist.append(btn)
                                    
                    
        label = Label(self.frame, text=label_text, fg=label_color, bg=label_bgcolor)
        label.grid(row=label_row, column=0, columnspan=self.columnspantextlabel, sticky=E)
        w.widgetslist.append(label)
        
        return w
        
    def fillHeader(self):
        import os, sys        
        self.master.title(self.programname)
        headertext = 'GUI for Xmipp %s \n Executed in directory: %s' % (self.programname, os.getcwd())
        self.l1 = Label(self.frame, text=headertext, fg=self.style.SectionTextColor, bg=self.style.LabelBgColor)
        self.l1.configure(wraplength=self.style.WrapLenght)
        self.l1.grid(row=self.getRow(), column=0, columnspan=6, sticky=E+W)
        self.citerow = self.getRow()        
        self.addSeparator(self.getRow())
            
    def browse(self, var, isFile):
        import tkFileDialog
        import os
        if isFile:
            filename = tkFileDialog.askopenfilename(title="Choose file")
        else:
            filename = tkFileDialog.askdirectory(title="Choose directory")
        if len(filename) > 0:
            var.set(os.path.relpath(filename))
         
    def scroll(self, event):
        if event.num == 5 or event.delta == -120:
            count = 1
        if event.num == 4 or event.delta == 120:
            count = -1
        self.canvas.yview("scroll", count,"units")
    #-------------------------------------------------------------------
    # Reading and parsing script
    #-------------------------------------------------------------------
    def readProtocolScript(self):
        begin_of_header = False
        end_of_header = False        
        f = open(self.scriptname, 'r')
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
        reVariable = re.compile('(\w+)\s*=\s*(True|False|".*"|\'.*\'|\d+|)')
        self.variablesDict = {}
        self.widgetslist = []
        lastSection = None
        
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
                    if index < count and self.header_lines[index].startswith('"""'):
                        while index < count and not helpStr.endswith('"""\n'):
                            line = self.header_lines[index]
                            helpStr += line
                            index += 1
                        v.help = helpStr
                        
                    if 'please_cite' in v.tags.keys():
                        self.citeslist.append(v.help.replace('"', ''))
                    else:
                        if index < count:
                            line = self.header_lines[index].strip()
                            match2 = reVariable.match(line)
                            if match2:
                                v.name, v.value = (match2.group(1).strip(), match2.group(2).strip())
                                #print "DEBUG_JM: v.name: '%s', v.value: '%s'" % (v.name, v.value)
                                self.variablesDict[v.name] = v
                                index += 1
                if is_section or v.name:
                    w = self.createWidget(v)
                    
                    
        #Update if citations found
        if len(self.citeslist):
            citetext = "If you publish results obtained with this protocol, please cite:\n"
            citetext += '\n'.join(self.citeslist)
            label = Label(self.frame, text=citetext, fg=self.style.CitationTextColor, bg=self.style.LabelBgColor)
            label.configure(wraplength=self.style.WrapLenght)
            label.grid(row=self.citerow, column=0, columnspan=5, sticky=EW)
                    
                

            
    #-------------------------------------------------------------------
    # GUI Events handling
    #-------------------------------------------------------------------           
        
    def close(self, event=""):
        self.master.destroy()
    
    def toggleExpertMode(self, event=""):
        self.expert_mode = not self.expert_mode
        
        if self.expert_mode:
            text = "Hide Expert Options"
        else:
            text = "Show Expert Options"
        self.btnExpert.config(text=text)
        self.checkVisibility()
    
    def save(self, event=""):
        print "* Saving script: %s" % self.scriptname
        f = open(self.scriptname, 'w')
        #f = sys.stdout
        f.writelines(self.pre_header_lines)
        for w in self.widgetslist:
            wlines = w.variable.getLines()
            f.writelines(wlines)
        f.writelines(sepLine + '# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE\n') 
        f.writelines(self.post_header_lines)
        f.close()
        os.chmod(self.scriptname, 0755)
    
    def confirmDeleteWorkingDir(self):
        if 'DoDeleteWorkingDir' in self.variablesDict and self.variablesDict['DoDeleteWorkingDir'].getValue() == 'True':
             return tkMessageBox.askyesno("Confirm DELETE", "Working dir '%s' will be DELETED. Do you want to continue?" % self.variablesDict['WorkingDir'].getValue())
        return True
    
    def validateProtocol(self):
        mod= loadModule(self.scriptname)
        #print dir(mod)
        if 'checkErrors' in dir(mod):
            errors = mod.checkErrors()
            if len(errors) > 0:
                tkMessageBox.showerror("Validation ERRORS",'\n'.join(errors))
                return False
        return True
    
    def saveExecute(self, event=""):
        self.save() 
        if self.confirmDeleteWorkingDir() and self.validateProtocol():
            print "EXECUTING"           
    
    def showHelp(self, helpmsg):
        import tkMessageBox
        tkMessageBox.showinfo("Help", helpmsg)
    
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
            self.style.Font.configure(size = new_size)
        self.resize()
        
    def checkVisibility(self, event=""):        
        for w in self.widgetslist:
            w.checkVisibility()    
    
    def fillButtons(self):
        row = self.getRow()
        self.addSeparator(row)
        row += 3
        self.addButton("Close", self.close, 0, row, 0, W)
        self.btnExpert = self.addButton("Show Expert Options", self.toggleExpertMode, 12, row, 1, EW)
        self.addButton("Save", self.save, 0, row, 3, W)
        self.addButton("Save & Execute", self.saveExecute, 7, row, 4, W)
        
    def addBindings(self):
        self.master.bind('<Alt_L><c>', self.close)
        self.master.bind('<Alt_L><o>', self.toggleExpertMode)
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
        
    def createGUI(self, script):
        self.init(script)        
        #self.master = Tk()
        #self.style = ProtocolStyle('config_gui')
        #self.style.createFonts()
        self.createBasicGUI()
        self.master.option_add("*Font", self.style.Font)
        self.columnspantextlabel = 3
        self.columntextentry = 3
        
        self.readProtocolScript()
        self.createScrollableCanvas()
        self.fillHeader()
        self.parseHeader()
        #self.fillWidgets()
                # Add bottom row buttons
    def fillGUI(self):
        self.fillButtons()
        self.addBindings()        
        self.checkVisibility()  
    
if __name__ == '__main__':
    script = sys.argv[1]  
    gui = ProtocolGUI()
    gui.createGUI(script)
    gui.launchGUI()
    
