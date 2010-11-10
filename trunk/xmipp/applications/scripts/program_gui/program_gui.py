#!/usr/bin/env python
"""/***************************************************************************
 *
 * Author: J.M. de la Rosa Trevin     jmdelarosa@cnb.csic.es
 *
 * Universidad Autonoma de Madrid
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
import sys, fileinput;
from Tkinter import *;  
import tkFont;
import tkMessageBox;
import tkFileDialog;

BORDER = 0;

class OptionWidget(Frame):
    def __init__(self, parent, optionName, default, subOptions=0):
        self.parent = parent;
        self.default = default;
        self.subparams = [];
        self.lastIndex = -1;
        Frame.__init__(self, parent, bd=BORDER, relief=SUNKEN);
        Label(self, text=str(optionName)).pack(side=LEFT, padx="2m");
        self.isFile = False;
        if subOptions == 0:
            entryWidth = max(6, len(default));
            iName = optionName.lower();
            self.isFile = iName.find("file") != -1 or iName.find("metadata") != -1;
            if self.isFile:
                entryWidth = 40;
            self.entryVar = StringVar();
            self.entryVar.set(str(default));
            self.entry = Entry(self,
              bg="white",
              width=entryWidth,
              highlightcolor="magenta",
              textvariable=self.entryVar,
              );
            self.entry.pack(side=LEFT, padx="1m");
            if self.isFile:
                self.browseButton = Button(self, text="Browse", command=self.browseFile);
                self.browseButton.pack(side=LEFT, padx="1m");
        else:
            self.menuTextVar = StringVar();
            self.maxWidth = 6;    
            self.menuSelectionVar = IntVar();
            self.menuSelectionVar.set(-1);
            self.menuButton = Menubutton(self,
                                         relief=RAISED,
                                         textvariable=self.menuTextVar);            
            self.menuButton.menu = Menu(self.menuButton);
            self.menuButton["menu"] = self.menuButton.menu;
            self.menuButton.pack(side=LEFT);
            
    def browseFile(self):
        filename = tkFileDialog.askopenfilename();
        if len(filename) > 0:
            self.entryVar.set(filename);
            
    def setState(self, optState):
        if len(self.subparams) == 0:
            self.entry.config(state=optState);
            if self.isFile:
                self.browseButton.config(state=optState);
        else:
            self.menuButton.config(state=optState);
            
    def addParam(self, param):
        self.subparams.append(param);
        self.maxWidth = max(self.maxWidth, len(param.paramName));
        self.menuButton.config(width=self.maxWidth);
        index = len(self.subparams) - 1;
        self.menuButton.menu.add_radiobutton(label=param.paramName,
                                             variable=self.menuSelectionVar,
                                             value=index,
                                             command=self.selectionChanged);
        if param.paramName != self.default:
            param.pack_forget();            
        else:
            self.menuTextVar.set(self.default);
            self.lastIndex = index;
        return Label(param, text="");        
    
    def selectionChanged(self):
        index = self.menuSelectionVar.get();
        self.menuTextVar.set(self.subparams[index].paramName);
        if self.lastIndex != -1:
            self.subparams[self.lastIndex].pack_forget();
        self.subparams[index].pack();
        self.lastIndex = index; 
             
class ParamWidget(Frame):
    def __init__(self, parent, paramName):
        self.comments = [];
        self.options = [];
        self.parent = parent;
        self.paramName = paramName;
        #The following flag will be used to pack the 
        #param in the same line or in the next one
        self.isLarge = False; 
        Frame.__init__(self, parent, bd=BORDER, relief=SUNKEN);
        self.pack(anchor="w", fill=X, side=TOP);
        self.normalFont = tkFont.Font(weight="normal");
        self.boldFont = tkFont.Font(weight="bold");  
        #Add label of radiobutton   
        self.label = parent.addParam(self);
        self.label.pack(side=LEFT);          
            
    def addOption(self, optionName, defaultValue, subOptions=0):
        self.lastOption = OptionWidget(self, optionName, defaultValue, subOptions)
        self.lastOption.pack(side=LEFT);
        self.options.append(self.lastOption);
        if self.lastOption.isFile:
            self.isLarge = True;
        return self.lastOption;
        
    def addSubParam(self, paramName):
        param = ParamWidget(self.lastOption, paramName);
        self.isLarge = True;
        return param;
        #self.lastOption.addParam("cosine");
    
    def endWithOptions(self):
        #Put a checkbutton if have no arguments           
        if len(self.options) == 0:
           self.checkButton = Checkbutton(self);
           self.checkButton.pack(side=LEFT);
        if len(self.comments) > 0:
            helpImage = PhotoImage(file="~/Download/help256.gif").subsample(12, 12);
            helpButton = Button(self, image=helpImage, command=self.buttonHelp_click);
            helpButton.image = helpImage;
            helpButton.pack(side=LEFT, padx="1m");  
        if not self.parent.single and len(self.parent.params) > 1:
            self.disable(); 
        if len(self.options) > 2:
            self.isLarge = True;            
         
        
    def addCommentLine(self, comment):
        self.comments.append(comment);
        
    def buttonHelp_click(self):
        msg = "\n".join(self.comments);
        tkMessageBox.showinfo(self.paramName + " help", msg);
        
    def enable(self):
        self.label.config(font=self.boldFont);
        if len(self.options) == 0:
            self.checkButton.config(state=NORMAL);
        else:
            for option in self.options:
                option.setState(NORMAL)
            
    def disable(self):
        self.label.config(font=self.normalFont);
        if len(self.options) == 0:
            self.checkButton.config(state=DISABLED);
        else:
            for option in self.options:
                option.setState(DISABLED);
                
    def getParamString(self):
        if len(self.options) == 0:
            pass;

class ParamsGroup(Frame):
    def __init__(self, parent, single):
        if single:
            border = 0;
        else:
            border = 1;
        Frame.__init__(self, parent, bd=border, relief=GROOVE);
        #self.pack(anchor=W, fill=X);        
        self.single = single;
        self.params = [];
        self.selectionVar = IntVar();
        self.selectionVar.set(0);
        self.lastSelected = 0;
        
    def addParam(self, param):
        self.params.append(param);
        if self.single:
            return Label(param, text=param.paramName, font=param.boldFont);
        else:
            return  Radiobutton(param, text=param.paramName,
                            variable=self.selectionVar,
                            value=len(self.params) - 1,
                            font=param.boldFont,
                            command=self.selectionChanged);
    
    def selectionChanged(self):
        self.params[self.lastSelected].disable();
        self.lastSelected = self.selectionVar.get();
        self.params[self.lastSelected].enable();
        
    def isLarge(self):
        return not self.single or self.params[0].isLarge;
    
    def getParamString(self):
        return self.params[self.lastSelected].paramName;
        
class SectionWidget(LabelFrame):
    def __init__(self, parent, sectionName):        
        self.sectionName = sectionName.strip();
        LabelFrame.__init__(self, parent, bd=1, relief=SUNKEN,
                            text=self.sectionName,
                            fg="#297739",
                            font=tkFont.Font(weight="bold"),
                            labelanchor="nw");
        #To check when to pack params in next line
        self.row = 0;
        self.col = 0;
        self.groups = [];
    def addGroup(self, group):
        #group.pack(fill=X, anchor=W);
        #return;
        #print "section: ", self.sectionName;
        #print "group: ", group.params[0].paramName#, " col: ", self.col, " row: ", self.row;
        #return;
        if group.isLarge():            
            group.grid(column=0, columnspan=2, row=self.row, 
                       sticky=NW, ipadx="2m");
            self.row = self.row + 1;
            self.col = 0;
        else:
            group.grid(column=self.col, row=self.row, 
                       sticky=NW, ipadx="2m");
            if self.col == 1:
                self.row = self.row + 1;
                self.col = 0;
            else:
                self.col = 1;
        self.groups.append(group);
    
class ProgramGUI(Frame):
    def __init__(self, parent):
        self.parent = parent;
        Frame.__init__(self, parent);        
        self.sections = [];
        # Create header
        self.createHeader();
        # Create params frame and add params
        self.params_frame = Frame(self);
        self.params_frame.pack(fill=X, padx="5m");
        self.createParams();
        self.createButtons();
        self.pack();
    
    def createHeader(self):
        '''Create header with program name and usage lines'''
                # Read the title from input 
        title = sys.stdin.readline();#readline();
        self.winfo_toplevel().title(title);
        self.progName = title.split("-")[1].strip();
        Label(self,
              text=self.progName,
              font=tkFont.Font(weight="bold", size="16"),
              fg="#292987", bd=BORDER, relief=SUNKEN
              ).pack(side=TOP, padx="5m"); 
        n = int(sys.stdin.readline());
        if n > 0:
            usageLabel = Label(self, justify=LEFT, bd=BORDER, relief=SUNKEN);
            usageText = "";
            for i in range(0, n):
                usageText += sys.stdin.readline();
            usageLabel.config(text=usageText);
            usageLabel.pack(side=TOP, padx="5m");
            
              
    def createParams(self):
        #self.createTestingSection();
        #return;
        '''Create parameters reading from standard input'''    
        # Params added from the standard input
        # this mecanism will be used for interprocess comunication
        # for each program could add theirs owns params widgets
        for line in fileinput.input():
            #print "line: ", line.strip();
            exec(line);
            
        
    def createTestingSection(self):
        # Create some params for testing
        section = self.addSection("Testing section");
        group = ParamsGroup(section, False);
        
        
        param = ParamWidget(group, "--input");
        param.endWithOptions();
        param = ParamWidget(group, "--ouput");
        param.addCommentLine("bla bla bla");
        param.addCommentLine("more bla bla bla ");
        param.addOption("a", "100", 0);
        param.addOption("b", "0", 0);
        param.endWithOptions();        
        param = ParamWidget(group, "--ouput222");
        param.addCommentLine("bla bla bla");
        param.addCommentLine("more bla bla bla ");
        param.addOption("a", "100", 0);
        param.addOption("b", "0.", 0);
      
        param.endWithOptions();
        section.addGroup(group)
        
        group = ParamsGroup(section, True);
        param = ParamWidget(group, "--ouput3333");
        param.addCommentLine("bla bla bla");
        param.addCommentLine("more bla bla bla ");
        param.addOption("a", "100", 0);
        #Example of option with subparams
        param.addOption("wave type", "sine", 2);
        subparam = param.addSubParam("cosine");
        subparam.addOption("x", "0", 0);
        subparam.addOption("y", "", 0);
        subparam = param.addSubParam("sine");
        subparam.addOption("z", "100");  
        param.endWithOptions();
        section.addGroup(group);
        
        group = ParamsGroup(section, True);
        param = ParamWidget(group, "--mask");
        param.addCommentLine("bla bla bla");
        param.endWithOptions();
        section.addGroup(group);
        
        group = ParamsGroup(section, True);
        param = ParamWidget(group, "--cookie");
        param.addOption("a", "100", 0);
        param.addCommentLine("bla bla bla");
        param.endWithOptions();
        section.addGroup(group);
        
    def createButtons(self):
        # Create buttons frame and add butttons
        self.buttons_frame = Frame(self);
        self.buttons_frame.pack(anchor="e", pady="3m");
        # Cancel button
        Button(self.buttons_frame,
               text="Cancel",
               width=6,
               command=self.buttonCancel_click
               ).pack(side=RIGHT, padx="1m", anchor="e");
        # Run button
        Button(self.buttons_frame,
               text="Run",
               width=6,
               command=self.buttonRun_click
               ).pack(side=RIGHT, padx="1m", anchor="e");
        Button(self.buttons_frame,
               text="Center",
               width=6,
               command=self.centerWindows
               ).pack(side=RIGHT, padx="1m", anchor="e");
    def addSection(self, sectionName):
        section = SectionWidget(self.params_frame, sectionName);
        self.sections.append(section);
        section.pack(padx="2m", side=TOP, fill=X);
        return section;
        #Label(self.params_frame, text=sectionName + "Test").pack();
        
    def buttonCancel_click(self):
        self.parent.destroy();
        
    def buttonRun_click(self):
        cmdStr = self.progName;
        for section in self.sections:
            for group in section.groups:
                cmdStr = cmdStr + " " + group.getParamString();
        print cmdStr;
        
    def centerWindows(self):
        root = self.parent
        w = root.winfo_screenwidth()
        h = root.winfo_screenheight()
        print w, h
        rootsize = tuple(int(_) for _ in root.geometry().split('+')[0].split('x'))
        print rootsize
        x = w/2 - rootsize[0]/2
        y = h/2 - rootsize[1]/2
        root.geometry("%dx%d+%d+%d" % (rootsize + (x, y)))
        
        
root = Tk();
#myapp = SampleApp(root)
program = ProgramGUI(root);
#program.centerWindows();

root.mainloop();

