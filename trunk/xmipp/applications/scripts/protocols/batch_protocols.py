#!/usr/bin/env python
'''
/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
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
 '''

import os
import shutil 
from Tkinter import *
import tkFont
from protlib_filesystem import getXmippPath
from protlib_base import *
from protlib_utils import getScriptPrefix
from protlib_gui import *
from protlib_gui_ext import *
import tkMessageBox

#Font
FontName = "Helvetica"
FontSize = 10

#TextColor
CitationTextColor = "dark olive green"
LabelTextColor = "black"
SectionTextColor = "blue4"

#Background Color
BgColor = "white"
LabelBgColor = BgColor
HighlightBgColor = BgColor
ButtonBgColor = "LightBlue"
ButtonActiveBgColor = "LightSkyBlue"
EntryBgColor = "lemon chiffon" 
ExpertLabelBgColor = "light salmon"

#Color
ListSelectColor = "DeepSkyBlue4"
BooleanSelectColor = "DeepSkyBlue4"

#Dimensions limits
MaxHeight = 800
MaxWidth = 800
MaxFontSize = 14
MinFontSize = 6


        
class XmippProjectGUI():  
    def __init__(self, project):
        self.project = project
        
    def deleteProject(self):
        if tkMessageBox.askyesno("DELETE confirmation", "You are going to DELETE all project data (runs, logs, results...Do you want to continue?"):
            self.project.clean()
            self.updateRunHistory("")#clean history and details
            
    def createMainMenu(self):
        self.menubar = Menu(self.root)
        self.fileMenu = Menu(self.root, tearoff=0)
        self.fileMenu.add_command(label="Delete project", command=self.deleteProject)
        self.fileMenu.add_command(label="Exit", command=self.onExit)
        self.menubar.add_cascade(label="File", menu=self.fileMenu)
        self.ToolbarButtonsDict = {}
        self.runButtonsDict = {}
        self.lastSelected = None
        self.lastRunSelected = None
        self.root.bind('<Configure>', self.unpostMenu)
        self.root.bind("<Unmap>", self.unpostMenu)
        self.root.bind("<Map>", self.unpostMenu)
   
    def addHeaderLabel(self, parent, text, row, col=0):
        '''Add a label to left toolbar'''
        label = Label(parent, text=text, font=self.LabelFont, fg=SectionTextColor)
        label.grid(row = row, column=col)
        return label
        
    def createToolbarButton(self, row, text, opts=[]):
        '''Add a button to left toolbar'''
        Font = tkFont.Font(family=FontName, size=FontSize-1, weight=tkFont.BOLD)
        btn = Button(self.toolbar, bd = 1, text=text, font=self.ButtonFont, relief=RAISED,
                         bg=ButtonBgColor, activebackground=ButtonBgColor)
        btn.grid(row = row, column = 0, sticky=W+E, pady=2, padx=5)
        list = ['a', 'b', 'c', 'd']
        i = 0
        if len(opts) > 0:
            menu = Menu(self.root, bg=ButtonBgColor, activebackground=ButtonBgColor, font=self.ButtonFont, tearoff=0)
            prots = [protDict.protocolDict[o] for o in opts]
            for p in prots:
                #Following is a bit tricky, its due Python notion of scope, a for does not define a new scope
                # and menu items command setting
                def item_command(prot): 
                    def new_command(): 
                        self.newProtocol(prot)
                    return new_command 
                menu.add_command(label=p.title, command=item_command(p))
            menu.bind("<Leave>", self.unpostMenu)
        self.ToolbarButtonsDict[text] = (btn, menu)
        btn.config(command=lambda:self.selectToolbarButton(text))

    def launchProtocolGUI(self, run, srcScript, saveCallback):
        top = Toplevel()
        gui = ProtocolGUI()
        gui.createGUI(srcScript, self.project, run, top, saveCallback)
        gui.fillGUI()
        gui.launchGUI()
        
    def loadProtocol(self, prot, srcProtAbsPath):
        protocol = prot.key
        protDir = getXmippPath('protocols')
        #suggest a new run_name        
        runName = self.project.projectDb.suggestRunName(protocol)
        dstAbsPath = self.project.getRunScriptFileName(protocol, runName)
        run = {
               'protocol_name':protocol, 
               'run_name': runName, 
               'script': dstAbsPath, 
               'comment': "my first run"
               }
        s, ss = getSectionByKey(prot)
        self.launchProtocolGUI(run, srcProtAbsPath, 
                               lambda: self.protocolSaveCallback(run, ss))
                
    def newProtocol(self, prot):
        protDir = getXmippPath('protocols')
        srcProtName = 'xmipp_protocol_%s.py' % prot.key
        srcProtDir = getXmippPath('protocols')
        srcProtAbsPath = os.path.join(protDir, srcProtName)
        self.loadProtocol(prot, srcProtAbsPath)
        
    def protocolSaveCallback(self, run, protGroups):
        if self.lastSelected in protGroups:
            self.updateRunHistory(self.lastSelected) 

    def updateRunHistory(self, protGroup): 
        self.runs = self.project.projectDb.selectRuns(protGroup)
        self.lbHist.delete(0, END)
        for run in self.runs:
            self.lbHist.insert(END, ('%s_%s' % (run['protocol_name'], run['run_name']),
                                run['last_modified']))   
        if len(self.runs) > 0:
                self.lbHist.selection_set(0)         

    #---------------- Functions related with Popup menu ----------------------   
    def lastPair(self):
        if self.lastSelected:
            return  self.ToolbarButtonsDict[self.lastSelected]
        return None
        
    def unpostMenu(self, event=None):
        if self.lastSelected:
            btn, menu = self.lastPair()
            menu.unpost()
            
    def postMenu(self, btn, menu):
        x, y, w = btn.winfo_x(), btn.winfo_y(), btn.winfo_width()
        xroot, yroot = self.root.winfo_x(), self.root.winfo_y()
        menu.post(xroot + x + w + 10, yroot + y)
        btn.config(bg=ButtonActiveBgColor, activebackground=ButtonActiveBgColor)
        
    def selectToolbarButton(self, text):
        btn, menu = self.ToolbarButtonsDict[text]

        if self.lastSelected and self.lastSelected != text:
            lastBtn, lastMenu = self.lastPair()
            lastBtn.config(bg=ButtonBgColor)
            lastMenu.unpost()            
        if self.lastSelected:
            self.postMenu(btn, menu)
        
        if self.lastSelected != text:
            self.project.config.set('project', 'lastselected', text)
            self.project.writeConfig()
            self.updateRunSelection(-1)
            self.updateRunHistory(text)            
        self.lastSelected = text  
            
    def updateRunSelection(self, index):
        state = NORMAL
        if index == -1:
            state = DISABLED
        for btn in self.runButtonsDict.values():
            btn.config(state=state)
            
        #def updateRunDetails(self, index):
        #    run = self.runs[index]
            
    def runSelectCallback(self, index):
        if index >= 0:
            self.lastRunSelected = self.runs[index]
        self.updateRunSelection(index)
            
        
    def runButtonClick(self, event=None):
        run = dict(zip(runColumns, self.lastRunSelected))
        s, ss = getSectionByKey(protDict.protocolDict[run['protocol_name']])
        if event == 'Edit':
            self.launchProtocolGUI(run, run['script'], 
                               lambda: self.protocolSaveCallback(run, ss))
        elif event == 'Copy':
            self.loadProtocol(protDict.protocolDict[run['protocol_name']], run['script'])
        elif event == "Delete":
            if tkMessageBox.askyesno("Confirm DELETE", "All data related to this run will be DELETED. Do you want to continue?"):
                self.project.deleteRun(run)
                self.updateRunHistory(self.lastSelected)
        

    def createToolbar(self):
        #Configure toolbar frame
        self.toolbar = Frame(self.frame, bd=2, relief=RIDGE)
        self.toolbar.grid(row=0, column=0, sticky=N+W+S, 
                          rowspan=5, padx=5, pady=5)
        #Create toolbar buttons
        i = 1
        for k, v in sections:
            self.addHeaderLabel(self.toolbar, k, i)
            i += 1
            for btn in v:
                self.createToolbarButton(i, btn[0], btn[1:])
                i += 1
                
    def addRunButton(self, frame, text, col, imageFilename=None):
        btnImage = None
        if imageFilename:
            try:
                imgPath = os.path.join(getXmippPath('resources'), imageFilename)
                btnImage = PhotoImage(file=imgPath)
            except TclError:
                pass
        
        if btnImage:
            btn = Button(frame, image=btnImage, bd=0, height=28, width=28)
            btn.image = btnImage
        else:
            btn = Button(frame, text=text, font=self.ButtonFont, bg=ButtonBgColor)
        btn.config(command=lambda:self.runButtonClick(text), 
                 activebackground=ButtonActiveBgColor)
        btn.grid(row=0, column=col)
        ToolTip(btn, text, 500)
        self.runButtonsDict[text] = btn
    
    def createRunHistory(self):
        label = self.addHeaderLabel(self.frame, 'History', 0, 1)
        #Button(self.frame, text="Edit").grid(row=0, column=2)
        frame = Frame(self.frame)
        frame.grid(row=0, column=2)
        self.addRunButton(frame, "Edit", 0, 'edit.gif')
        self.addRunButton(frame, "Copy", 1, 'copy.gif')
        self.addRunButton(frame, "Visualize", 2, 'visualize.gif')
        self.addRunButton(frame, "Delete", 3, 'delete.gif')
        self.frameHist = Frame(self.frame)
        self.frameHist.grid(row=2, column=1, sticky=NSEW, columnspan=2)
        self.lbHist = MultiListbox(self.frameHist, (('Run', 40), ('Modified', 20)))
        self.lbHist.SelectCallback = self.runSelectCallback
        self.lbHist.DoubleClickCallback = lambda:self.runButtonClick("Edit")
        self.lbHist.AllowSort = False        
        #self.lbHist.pack()
        
    def createRunDetails(self):
        #Create RUN details
        self.addHeaderLabel(self.frame, 'Details', 3, 1)
        self.frameDetails = Frame(self.frame, bg=BgColor, bd=1, relief=RIDGE)
        self.frameDetails.grid(row=4, column=1,sticky=NSEW, columnspan=2)

    def createGUI(self, root=None):
        if not root:
            root = Tk()
        self.root = root
        root.withdraw() # Hide the windows for centering
        self.root.title("Xmipp Protocols")
        self.createMainMenu()
        #Create a main frame that contains all other widgets
        self.frame = Frame(self.root)
        self.frame.pack(fill=BOTH)
        self.frame.columnconfigure(0, minsize=150, weight=1)
        self.frame.columnconfigure(1, minsize=300, weight=2)
        #self.frame.columnconfigure(2, minsize=300, weight=2)
        self.frame.rowconfigure(2, minsize=50, weight=1)
        self.frame.rowconfigure(4, minsize=50, weight=1)
        
        # Create some fonts for later use
        self.ButtonFont = tkFont.Font(family=FontName, size=FontSize, weight=tkFont.BOLD)
        self.LabelFont = tkFont.Font(family=FontName, size=FontSize+1, weight=tkFont.BOLD)
        
        self.createToolbar()
        self.createRunHistory()
        self.createRunDetails()

        self.root.config(menu=self.menubar)
        #select lastSelected
        if self.project.config.has_option('project', 'lastselected'):
            self.selectToolbarButton(self.project.config.get('project', 'lastselected'))
    
    def launchGUI(self, center=True):
        if center:
            self.root.update_idletasks()
            centerWindows(self.root)
        self.root.deiconify()
        self.root.mainloop()
       
    def onExit(self):
        self.root.destroy()


if __name__ == '__main__':
    import sys
    dir = os.getcwd()
    project = XmippProject(dir)    
    
    if len(sys.argv) > 1:
        # Launch a protocol directly
        from protocol_gui import *
        script = sys.argv[1]
        project.load()  
        gui = ProtocolGUI()
        gui.createGUI(script)
        gui.launchGUI()
     
    else: #lauch project     
        projectCfg = '.project.cfg'
        if not os.path.exists(projectCfg):
            print 'You are in directory: ', dir
            answer = raw_input('Do you want to create a new xmipp_protocols PROJECT in this folder? [Y/n]:')
            if not answer or answer.lower() == 'y':
                project.create()
        else:
            project.load()
        gui = XmippProjectGUI(project)
        gui.createGUI()
        gui.launchGUI()
    
