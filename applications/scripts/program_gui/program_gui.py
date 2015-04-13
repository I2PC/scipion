#!/usr/bin/env python
"""/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
 """
import sys, os, fileinput, re;
from Tkinter import *;  
import tkFont;
import tkMessageBox;
import tkFileDialog;
from subprocess import call
from protlib_filesystem import getXmippPath

BORDER = 0;
PROGRAM_BLUE = "#292987";
SECTION_GREEN = "#297739";
LARGE_FONT_SIZE = "16"
NORMAL_FONT_SIZE = "12"
SMALL_FONT_SIZE = "10" 

def isEntryFile(name):
    return name.find("file") != -1 or name.find("metadata") != -1 or name.find("directory") != -1

def getEntryWith(name, default=None):
    if isEntryFile(name):
        return 40;
    if name.endswith("s") or name.find("expression") != -1: # get a plural optionName
        return 20;
    if default:
        return max(10, len(default));
    return 10;

def createImageButton(parent, text, image, cmd):
    try:
        helpImage = PhotoImage(file = getXmippPath(os.path.join("resources", image)))
        btn = Button(parent, image=helpImage, command=cmd, bd=0)
        btn.image = helpImage;
    except TclError:
        btn = Button(parent, text=text, bg="#a7dce1", activebackground="#72b9bf",
                            command=cmd)
    return btn
        
class OptionWidget(Frame):
    def __init__(self, parent, optionName, default, subOptions=0):
        self.parent = parent;
        self.default = default;
        self.subparams = [];
        self.lastIndex = -1;
        Frame.__init__(self, parent, bd=BORDER, relief=SUNKEN)
        Label(self, text=str(optionName)).pack(side=LEFT, padx="2m")
        self.isFile = False;
        if subOptions == 0:
            iName = optionName.lower()
            self.isFile = isEntryFile(iName)
            entryWidth = getEntryWith(iName, default)
            self.entryVar = StringVar()
            self.entryVar.set(str(default))
            self.entry = Entry(self,
              bg="white",
              width=entryWidth,
              highlightcolor=PROGRAM_BLUE,
              textvariable=self.entryVar,
              )
            self.entry.pack(side=LEFT, padx="1m")
            if self.isFile:
                self.browseButton = createImageButton(self, "Browse", "fileopen.gif", self.browseFile)
                self.browseButton.pack(side=LEFT, padx="1m")
        else:
            self.menuTextVar = StringVar()
            self.maxWidth = 6;    
            self.menuSelectionVar = IntVar()
            self.menuSelectionVar.set(-1)
            self.menuButton = Menubutton(self,
                                         relief=RAISED,
                                         textvariable=self.menuTextVar)            
            self.menuButton.menu = Menu(self.menuButton)
            self.menuButton["menu"] = self.menuButton.menu;
            self.menuButton.pack(side=LEFT)
            
    def browseFile(self):
        filename = tkFileDialog.askopenfilename()
        if len(filename) > 0:
            cwd = os.getcwd() + os.sep;
            if filename.startswith(cwd):
                filename = filename.replace(cwd, "")
            self.entryVar.set(filename)
            
    def setState(self, optState):
        if len(self.subparams) == 0:
            self.entry.config(state=optState)
            if self.isFile:
                self.browseButton.config(state=optState)
        else:
            self.menuButton.config(state=optState)
            
    def addParam(self, param):
        param.isSubparam = True;
        self.subparams.append(param)
        self.maxWidth = max(self.maxWidth, len(param.paramName))
        self.menuButton.config(width=self.maxWidth)
        index = len(self.subparams) - 1;
        self.menuButton.menu.add_radiobutton(label=param.paramName,
                                             variable=self.menuSelectionVar,
                                             value=index,
                                             command=self.selectionChanged)
        if param.paramName != self.default:
            param.pack_forget()            
        else:
            self.menuTextVar.set(self.default)
            self.lastIndex = index;
        return Label(param, text="")        
    
    def selectionChanged(self):
        index = self.menuSelectionVar.get()
        self.menuTextVar.set(self.subparams[index].paramName)
        if self.lastIndex != -1:
            self.subparams[self.lastIndex].pack_forget()
        self.subparams[index].pack()
        self.lastIndex = index; 
        
    def getCmdString(self, showDefault):
        if len(self.subparams) == 0:
            if showDefault or not self.isDefault():
                return self.entryVar.get().strip()
            return "";
        else:
            if showDefault or not self.isDefault():
                return self.subparams[self.lastIndex].getCmdString(showDefault)
            return "";
    
    def isDefault(self):
        if len(self.subparams) == 0:
            return self.default == self.entryVar.get().strip()
        else:
            if self.lastIndex == -1:
                return True;                
            if self.default != self.subparams[self.lastIndex].paramName:
                return False;
            return self.subparams[self.lastIndex].isDefault()
        
    def getHelp(self):
        if len(self.subparams) == 0:
            return "";
        return self.subparams[self.lastIndex].getHelp()
        
             
class ParamWidget(Frame):
    def __init__(self, parent, paramName):
        self.comments = [];
        self.options = [];
        self.parent = parent;
        self.paramName = paramName;
        self.isSubparam = False;
        self.notOptional = False;
        #The following flag will be used to pack the 
        #param in the same line or in the next one
        self.isLarge = False; 
        Frame.__init__(self, parent, bd=BORDER, relief=SUNKEN)
        self.pack(anchor="w", fill=X, side=TOP)
        self.normalFont = tkFont.Font(weight="normal", size=NORMAL_FONT_SIZE)
        self.boldFont = tkFont.Font(weight="bold", size=NORMAL_FONT_SIZE)  
        #Add label of radiobutton   
        self.label = parent.addParam(self)
        self.label.pack(side=LEFT) 
            
    def addOption(self, optionName, defaultValue, subOptions=0):
        self.lastOption = OptionWidget(self, optionName, defaultValue, subOptions)
        self.lastOption.pack(side=LEFT)
        self.options.append(self.lastOption)
        if self.lastOption.isFile:
            self.isLarge = True;
        return self.lastOption;
        
    def addSubParam(self, paramName):
        param = ParamWidget(self.lastOption, paramName)
        self.isLarge = True;
        return param;
        #self.lastOption.addParam("cosine")
    
    def endWithOptions(self):
        #Put a checkbutton if have no arguments           
        if len(self.options) == 0 and self.parent.single:
            self.checkVar = IntVar()
            self.checkButton = Checkbutton(self, variable=self.checkVar)
            self.checkButton.pack(side=LEFT)
        if len(self.comments) > 0:
            helpButton = createImageButton(self, "Help", "help.gif", self.buttonHelp_click)
            helpButton.pack(side=LEFT, padx="1m")              
        if not self.parent.single and len(self.parent.params) > 1:
            self.disable() 
        if len(self.options) > 2:
            self.isLarge = True;
        if self.notOptional:
            self.label.config(fg="#931212",
                              highlightcolor="#931212",
                              activeforeground="#931212")
    
    def addCommentLine(self, comment):
        self.comments.append(comment)
        
    def buttonHelp_click(self):
        msg = self.getHelp() + "\n";
        for option in self.options:
            msg = msg + "\n" + option.getHelp()
        #print msg;
        showInfo(self.paramName + " help", msg)
        
    def enable(self):
        self.label.config(font=self.boldFont)
        if len(self.options) == 0 and self.parent.single:
            self.checkButton.config(state=NORMAL)
        else:
            for option in self.options:
                option.setState(NORMAL)
            
    def disable(self):
        self.label.config(font=self.normalFont)
        if len(self.options) == 0 and self.parent.single:
            self.checkButton.config(state=DISABLED)
        else:
            for option in self.options:
                option.setState(DISABLED)
                
    def getHelp(self):
        return "\n".join(self.comments)
    
    def getCmdString(self, showDefault):
        if len(self.options) == 0:
            if self.isSubparam or not self.parent.single or self.checkVar.get() == 1:
                return " " + self.paramName;
            return "";
        else:
            optStr = "";
            if showDefault or not self.isDefault():                
                for option in reversed(self.options):
                    if not option.isDefault():
                        showDefault = True;
                    optStr = option.getCmdString(showDefault) + " " + optStr;
                optStr = self.paramName + " " +  optStr;
            return optStr;
            
    def isDefault(self):
        for option in self.options:
            if not option.isDefault():
                return False;
        return True;
    
    def check(self):
        if len(self.options) == 0:
            return "";
        if self.notOptional and self.isDefault():
            return "You should provide parameter: " + self.paramName + "\n";
        return "";
        

class ParamsGroup(Frame):
    def __init__(self, parent, single):
        if single:
            border = 0;
        else:
            border = 1;
        Frame.__init__(self, parent, bd=border, relief=GROOVE)
        #self.pack(anchor=W, fill=X)        
        self.single = single;
        self.params = [];
        self.selectionVar = IntVar()
        self.selectionVar.set(0)
        self.lastSelected = 0;
        
    def addParam(self, param):
        self.params.append(param)
        if self.single:
            return Label(param, text=param.paramName, font=param.boldFont)
        else:
            return  Radiobutton(param, text=param.paramName,
                            variable=self.selectionVar,
                            value=len(self.params) - 1,
                            font=param.boldFont,
                            command=self.selectionChanged)
    
    def selectionChanged(self):
        self.params[self.lastSelected].disable()
        self.lastSelected = self.selectionVar.get()
        self.params[self.lastSelected].enable()
        
    def isLarge(self):
        return not self.single or self.params[0].isLarge;
    
    def getCmdString(self, showDefault = False):
        return self.params[self.lastSelected].getCmdString(showDefault)
    
    def check(self):
        return self.params[self.lastSelected].check()
        
class SectionWidget(LabelFrame):
    def __init__(self, parent, sectionName):        
        self.sectionName = sectionName.strip()
        LabelFrame.__init__(self, parent, bd=1, relief=SUNKEN,
                            text=self.sectionName,
                            fg=SECTION_GREEN,
                            font=tkFont.Font(weight="bold", size=NORMAL_FONT_SIZE),
                            labelanchor="nw")
        #To check when to pack params in next line
        self.row = 0;
        self.col = 0;
        self.groups = [];
    def addGroup(self, group):
        self.groups.append(group)
        #group.pack(fill=X, anchor=W)
        #return;
        if group.isLarge(): 
            if self.col == 1:
                self.row = self.row + 1;           
            group.grid(column=0, columnspan=2, row=self.row,
                       sticky=NW, ipadx="2m")
            self.row = self.row + 1;
            self.col = 0;
        else:
            group.grid(column=self.col, row=self.row,
                       sticky=NW, ipadx="2m")
            if self.col == 1:
                self.row = self.row + 1;
                self.col = 0;
            else:
                self.col = 1;
        #print "section: ", self.sectionName;
        #print "group: ", group.params[0].paramName, " col: ", self.col, " row: ", self.row;
        #print "large: ", group.isLarge()
        
    
class ProgramGUI(Frame):
    def __init__(self, parent):
        self.parent = parent;
        Frame.__init__(self, parent)        
        self.sections = [];
        # Create header
        self.createHeader()
        # Create params frame and add params
        self.params_frame = Frame(self)
        self.params_frame.pack(fill=X, padx="5m")
        self.createParams()
        self.createButtons()   
        self.pack()     
    
    def createHeader(self):
        '''Create header with program name and usage lines'''
        # Read the title from input 
        title = sys.stdin.readline()#readline()
        self.winfo_toplevel().title(title)
        self.progName = title.split("-")[1].strip()
        Label(self,
              text=self.progName,
              font=tkFont.Font(weight="bold", size=LARGE_FONT_SIZE),
              fg=PROGRAM_BLUE, bd=BORDER, relief=SUNKEN
              ).pack(side=TOP, padx="5m") 
        n = int(sys.stdin.readline())
        if n > 0:
            usageLabel = Label(self, justify=LEFT, anchor=W, wraplength=500, bd=BORDER, relief=SUNKEN)
            usageText = "";
            for i in range(0, n):
                usageLine = sys.stdin.readline()
                usageText += usageLine.strip()     
            usageLabel.config(text=usageText)
            usageLabel.pack(side=TOP, padx="5m")
            
              
    def createParams(self):
        #self.createTestingSection()
        #return;
        '''Create parameters reading from standard input'''    
        # Params added from the standard input
        # this mecanism will be used for interprocess comunication
        # for each program could add theirs owns params widgets
        for line in fileinput.input():
            #print "line: ", line.strip()
            exec(line)
            
        
    def createTestingSection(self):
        # Create some params for testing
        section = self.addSection("Testing section")
        group = ParamsGroup(section, False)
        param = ParamWidget(group, "--input")
        param.endWithOptions()
        param = ParamWidget(group, "--ouput")
        param.addCommentLine("bla bla bla")
        param.addCommentLine("more bla bla bla ")
        param.addOption("a", "100", 0)
        param.addOption("b", "0", 0)
        param.endWithOptions()        
        param = ParamWidget(group, "--ouput222")
        param.addCommentLine("bla bla bla")
        param.addCommentLine("more bla bla bla ")
        param.addOption("a", "100", 0)
        param.addOption("b", "0.", 0)
      
        param.endWithOptions()
        section.addGroup(group)
        
        group = ParamsGroup(section, True)
        param = ParamWidget(group, "--ouput3333")
        param.addCommentLine("bla bla bla")
        param.addCommentLine("more bla bla bla ")
        param.addOption("a", "100", 0)
        #Example of option with subparams
        param.addOption("wave type", "sine", 2)
        subparam = param.addSubParam("cosine")
        subparam.addOption("x", "0", 0)
        subparam.addOption("y", "", 0)
        subparam = param.addSubParam("sine")
        subparam.addOption("z", "100")  
        param.endWithOptions()
        section.addGroup(group)
        
        group = ParamsGroup(section, True)
        param = ParamWidget(group, "--mask")
        param.addCommentLine("bla bla bla")
        param.endWithOptions()
        section.addGroup(group)
        
        group = ParamsGroup(section, True)
        param = ParamWidget(group, "--cookie")
        param.addOption("a", "100", 0)
        param.addCommentLine("bla bla bla")
        param.endWithOptions()
        section.addGroup(group)
        
    def createButtons(self):
        # Create buttons frame and add butttons
        self.buttons_frame = Frame(self)
        self.buttons_frame.pack(anchor="e", pady="3m", fill=X)
        # Cancel button
        Button(self.buttons_frame,
               text="Cancel",
               width=6,
               command=self.buttonCancel_click
               ).pack(side=RIGHT, padx="1m", anchor="e")
        # Run button
        Button(self.buttons_frame,
               text="Run",
               width=6,
               command=self.buttonRun_click
               ).pack(side=RIGHT, padx="1m", anchor="e")
               
        Button(self.buttons_frame,
               text="Show command",
               width=10,
               command=self.buttonShowCmd_click
               ).pack(side=RIGHT, padx="1m", anchor="e")


    def addSection(self, sectionName):
        section = SectionWidget(self.params_frame, sectionName)
        self.sections.append(section)
        section.pack(padx="2m", side=TOP, fill=X)
        return section;
        
    def check(self):
        errors = "";
        for section in self.sections:
            for group in section.groups:
                errors += group.check()
        if errors == "":
            return True;
        tkMessageBox.showerror(self.progName + " errors", errors)
        return False;
    
    def getCmd(self):
        cmdStr = self.progName + " ";
        for section in self.sections:
            for group in section.groups:
                cmdStr = cmdStr + group.getCmdString()
        return cmdStr;
        
    def buttonCancel_click(self):
        self.parent.destroy()
        
    def buttonRun_click(self):
        if self.check():
            cmdStr = self.getCmd()
            self.parent.destroy()
            print cmdStr;
            os.system(cmdStr)
            
    def copyCmdToClipboard(self):
        cmdStr = self.getCmd()
        self.clipboard_clear()
        self.clipboard_append(cmdStr)
        
    def buttonShowCmd_click(self):
        if self.check():
            cmdStr = self.getCmd()
            cmdWin = Toplevel(self.parent)
            cmdWin.title(self.progName + " command")
            text = Text(cmdWin, width=70, height=7, background="white")
            text.grid(row=0, column=0, columnspan=3)
            text.insert('1.0', cmdStr)
            Button(cmdWin, text="Copy to clipboard", command=self.copyCmdToClipboard).grid(row=1, column=1)
            Button(cmdWin, text="Close", width=6, command=cmdWin.destroy).grid(row=1, column=2)
            # Calculate window size to center command window
            pw, ph, px, py = getGeometry(self.parent)
            cmdWin.update_idletasks()
            centerWindows(cmdWin, pw, ph, px, py)
             
def getGeometry(win):
    return win.winfo_reqwidth(), win.winfo_reqheight(), win.winfo_x(), win.winfo_y()

def centerWindows(win, w, h, x, y):
    gw, gh, gx, gy = getGeometry(win)
    x = x + w / 2 - gw / 2
    y = y + h / 2 - gh / 2
    win.geometry("%dx%d+%d+%d" % (gw, gh, x, y))
        
        
root = Tk()
root.withdraw()
w = root.winfo_screenwidth()
h = root.winfo_screenheight()
program = ProgramGUI(root)
root.update_idletasks();
centerWindows(root, w, h, 0, 0)
root.deiconify()
root.mainloop()


