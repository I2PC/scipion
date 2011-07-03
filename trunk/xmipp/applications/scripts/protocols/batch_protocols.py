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
from config import *
from protlib_gui import *
from protlib_gui_ext import MultiListbox
  

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
        
    def createMainMenu(self):
        self.menubar = Menu(self.root)
        self.fileMenu = Menu(self.root, tearoff=0)
        self.fileMenu.add_command(label="Exit", command=self.onExit)
        self.menubar.add_cascade(label="File", menu=self.fileMenu)
        self.btnFrameDict = {}
   
    def addHeaderLabel(self, parent, text, row, col=0):
        '''Add a label to left toolbar'''
        label = Label(parent, text=text, font=self.LabelFont, fg=SectionTextColor)
        label.grid(row = row, column=col)
        
    def addLaunchButton(self, o, btnFrame, row, Font):
        label = Label(btnFrame, text=o, font=Font)
        label.grid(row=2*row, column=0, sticky=W, padx=5)
        btnLaunch = Button(btnFrame, text='New', font=Font, relief=RAISED,
                         bg=ButtonBgColor, activebackground=ButtonBgColor, command=lambda:self.newProtocol(o))
        btnLaunch.grid(row=2*row+1, column=0, padx=5, pady=5, sticky=E)
        
    def addTbButton(self, row, text, opts=[]):
        '''Add a button to left toolbar'''
        Font = tkFont.Font(family=FontName, size=FontSize-1, weight=tkFont.BOLD)
        btn = Button(self.toolbar, bd = 1, text=text, font=self.ButtonFont, relief=RAISED,
                         bg=ButtonBgColor, activebackground=ButtonBgColor)
        btn.grid(row = row, column = 0, sticky=W+E, pady=2, padx=5)
        
        btnFrame = Frame(self.frame, width=100, bd=1, relief=GROOVE)
        btnFrame.columnconfigure(0, minsize=150)
        label = Label(btnFrame, text='Protocols:', fg=SectionTextColor, font=Font)
        label.grid(row=0, column=0, pady=5)
        if len(opts) > 0:
            i = 0
            for o in opts:
                i += 1
                self.addLaunchButton(o, btnFrame, i, Font)
        else:
            self.addLaunchButton(text, btnFrame, 1, Font)
               
        self.btnFrameDict[text] = btnFrame 
        btn.config(command=lambda:self.menuPick(text))

    def newProtocol(self, protKey):
        protocol = launchDict[protKey]
        protDir = getXmippPath('protocols')
        srcProtName = 'xmipp_protocol_%s.py' % protocol
        srcProtDir = getXmippPath('protocols')
        srcProtAbsPath = os.path.join(protDir, srcProtName)

        #suggest a new run_name        
        lastRunName = self.project.projectDb.getLastRunName(protocol)        
        prefix, suffix  = getScriptPrefix(lastRunName)
        n = 1
        if suffix:
            n = int(suffix) + 1
        runName = "%s_%03d" % (prefix, n)
        dstAbsPath = os.path.join(self.project.runsDir, 'xmipp_protocol_%s_%s.py' % (protocol, runName))
#        print "Copying %s to %s" % (srcProtAbsPath, dstAbsPath)
#        shutil.copy(srcProtAbsPath, dstAbsPath)
        run = {
               'protocol_name':protocol, 
               'run_name': runName, 
               'script': dstAbsPath, 
               'comment': "my first run"
               }
        
        
        top = Toplevel()
        gui = ProtocolGUI()
        gui.createGUI(srcProtAbsPath, dstAbsPath, top)
        gui.fillGUI()
        #set the suggested runName
        gui.variablesDict['RunName'].setValue(runName)
        s, ss = getSection(protKey)
        gui.saveCallback = lambda: self.protocolSaveCallback(run, ss)
        gui.launchGUI()
        #os.system('python %s %s &' % (os.path.join(protDir, 'xmipp_protocol_gui.py'), dstAbsPath))
        
    def protocolSaveCallback(self, run, protGroup):
        self.project.projectDb.insertRun(run)
        
        if self.btnFrameDict[protGroup] == self.lastSelectedFrame:
            self.updateRunHistory(protGroup) 

    def updateRunHistory(self, protGroup): 
        self.runs = self.project.projectDb.selectRuns(protGroup)
        self.lbHist.delete(0, END)
        for run in self.runs:
            self.lbHist.insert(END, ('%s_%s' % (run['protocol_name'], run['run_name']),
                                run['last_modified']))
    def updateRunDetails(self, index):
        run = self.runs[index]
        
                
    
    def menuPick(self, text):
        frame = self.btnFrameDict[text]
        if self.lastSelectedFrame and frame != self.lastSelectedFrame:
            self.lastSelectedFrame.grid_remove()
        if frame != self.lastSelectedFrame:
            self.lastSelectedFrame = frame
            #self.lastSelectedFrame.grid(row=0, column=1, sticky=W+E+N, padx=10, pady=10, rowspan=4)
            self.project.config.set('project', 'lastselected', text)
            self.project.writeConfig()
            self.updateRunHistory(text)
            
    def runSelectCallback(self, index):
        print self.runs[index]['run_name']

    def createGUI(self, root=None):
        if not root:
            root = Tk()
        self.root = root
        self.root.title("Xmipp Protocols")
        self.createMainMenu()
        self.lastSelectedFrame = None
        #Create a main frame that contains all other widgets
        self.frame = Frame(self.root)
        self.frame.pack(fill=BOTH)
        self.frame.columnconfigure(0, minsize=150, weight=1)
        self.frame.columnconfigure(1, minsize=300, weight=2)
        #self.frame.columnconfigure(2, minsize=300, weight=2)
        self.frame.rowconfigure(1, minsize=50, weight=1)
        self.frame.rowconfigure(3, minsize=50, weight=1)
        
        #Configure toolbar frame
        self.toolbar = Frame(self.frame, bd=2, relief=RIDGE)
        self.toolbar.grid(row=0, column=0, sticky=N+W+S, 
                          rowspan=4, padx=5, pady=5)
        # Create some fonts for later use
        self.ButtonFont = tkFont.Font(family=FontName, size=FontSize, weight=tkFont.BOLD)
        self.LabelFont = tkFont.Font(family=FontName, size=FontSize+1, weight=tkFont.BOLD)
        #Create toolbar buttons
        i = 1
        for k, v in sections:
            self.addHeaderLabel(self.toolbar, k, i)
            i += 1
            for btn in v:
                self.addTbButton(i, btn[0], btn[1:])
                i += 1
        #Create RUNS history 
        self.addHeaderLabel(self.frame, 'History', 0, 1)
        self.frameHist = Frame(self.frame)
        self.frameHist.grid(row=1, column=1, sticky=N+W+E+S)
        self.lbHist = MultiListbox(self.frameHist, (('Run', 40), ('Modified', 20)))
        self.lbHist.SelectCallback = self.runSelectCallback
        self.lbHist.AllowSort = False
        self.lbHist.pack()
        #Create RUN details
        self.addHeaderLabel(self.frame, 'Details', 2, 1)
        self.frameDetails = Text(self.frame, bg=BgColor, height=10)
        self.frameDetails.grid(row=3, column=1)
        #self.detailsFrame = Frame(canvas)
        #self.detailsFrame.grid(row=0, column=0)
        
        self.root.config(menu=self.menubar)
        #select lastSelected
        if self.project.config.has_option('project', 'lastselected'):
            self.menuPick(self.project.config.get('project', 'lastselected'))
    
    def launchGUI(self):
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
    