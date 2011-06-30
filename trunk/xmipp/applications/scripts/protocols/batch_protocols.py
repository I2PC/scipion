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
   
    def addTbLabel(self, text, row):
        '''Add a label to left toolbar'''
        Font = tkFont.Font(family=FontName, size=FontSize+1, weight=tkFont.BOLD)
        label = Label(self.toolbar, text=text, font=Font, fg=SectionTextColor)
        label.grid(row = row, column=0)
        
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

    def newProtocol(self, btnName):
        protocol = launchDict[btnName]
        protDir = getXmippPath('protocols')
        srcProtName = 'xmipp_protocol_%s.py' % protocol
        srcProtDir = getXmippPath('protocols')
        srcProtAbsPath = os.path.join(protDir, srcProtName)
        
        lastRunName = self.project.projectDb.getLastRunName(protocol)
        
        g = getScriptPrefix(lastRunName)
        print g
        prefix, suffix =g
        n = 1
        if suffix:
            n = int(suffix) + 1
        runName = "%s_%03d" % (prefix, n)
        dstAbsPath = os.path.join(self.project.runsDir, 'xmipp_protocol_%s_%s.py' % (protocol, runName))
        print "Copying %s to %s" % (srcProtAbsPath, dstAbsPath)
        shutil.copy(srcProtAbsPath, dstAbsPath)
        run = {
               'protocol_name':protocol, 
               'run_name': runName, 
               'script': dstAbsPath, 
               'comment': "my first run"
               }
        self.project.projectDb.insertRun(run)
        
        top = TopLevel()
        gui = ProtocolGUI()
        gui.createGUI(script, top)
        gui.fillGUI()
        gui.launchGUI()
        #os.system('python %s %s &' % (os.path.join(protDir, 'xmipp_protocol_gui.py'), dstAbsPath))
        
    def launchProtocol(self, btnName):
        protName = 'xmipp_protocol_%s' % launchDict[btnName]
        protDestName = os.path.join(self.project.runsDir, protName)
        protDir = getXmippPath('protocols')
        
        if not os.path.exists(protDestName):
            protAbsPath = os.path.join(protDir, protName)
            shutil.copy(protAbsPath, protDestName)
        run = XmippProjectRun(protName, protDestName, "my first run")
        self.project.projectDb.insertRun(run)
        os.system('python %s %s &' % (os.path.join(protDir, 'xmipp_protocol_gui.py'), protDestName))
        
    def updateHistory(self, runs): 
        #TODO: this can be done in a better way
        for w in self.histFrame.grid_slaves():
            w.destroy()
        
        label = Label(self.histFrame, text="History", bg=BgColor, fg=SectionTextColor, font=self.ButtonFont)
        label.grid(row=0, column=0, columnspan=2)
        row = 1
        for run in runs:
            label = Label(self.histFrame, text=str(run['last_modified']), bg=BgColor)
            label.grid(row=row, column=1, sticky=W)
            label = Label(self.histFrame, text="%s_%s" % (run['protocol_name'], run['run_name']), 
                          bg="white", font=self.ButtonFont)
            label.grid(row=row, column=0, sticky=W, padx=15)
            row += 1       
        
    def menuPick(self, text):
        frame = self.btnFrameDict[text]
        if self.lastSelectedFrame and frame != self.lastSelectedFrame:
            self.lastSelectedFrame.grid_remove()
        if frame != self.lastSelectedFrame:
            self.lastSelectedFrame = frame
            self.lastSelectedFrame.grid(row=0, column=1, sticky=W+E+N, padx=10, pady=10)
            self.project.config.set('project', 'lastselected', text)
            self.project.writeConfig()
            runs = self.project.projectDb.selectRuns(text)
            self.updateHistory(runs)
            

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
        self.frame.columnconfigure(1, minsize=170, weight=1)
        self.frame.columnconfigure(2, minsize=300, weight=2)
        self.frame.rowconfigure(0, minsize=150)
        self.frame.rowconfigure(1, minsize=50)
        
        #Configure toolbar frame
        self.toolbar = Frame(self.frame, bd=2, relief=RIDGE)
        self.toolbar.grid(row=0, column=0, sticky=N+W+S, rowspan=2, padx=5, pady=5)#side=LEFT, fill=Y)
        # Create buttons
        self.ButtonFont = tkFont.Font(family=FontName, size=FontSize, weight=tkFont.BOLD)
        i = 1
        for k, v in sections:
            self.addTbLabel(k, i)
            i += 1
            for btn in v:
                self.addTbButton(i, btn[0], btn[1:])
                i += 1
            
        
        vscrollbar = Scrollbar(self.frame)
        vscrollbar.grid(row=0, column=3, sticky=N + S)
        canvas = Canvas(self.frame, width=50, height=150, bg=BgColor, bd=2,
                        yscrollcommand=vscrollbar.set, relief=RIDGE)
        canvas.grid(row=0, column=2, padx=5, pady=5, sticky=N+W+E+S)
        self.histFrame = Frame(canvas, bg=BgColor)
        self.histFrame.grid(row=0, column=0, pady=10, padx=10, sticky=N+W+E+S)
        canvas = Canvas(self.frame, height=50, bg=BgColor, bd=2, relief=RIDGE)
        canvas.grid(row=1, column=1, columnspan=3, sticky=S+W+E+N, padx=5, pady=5)
        self.detailsFrame = Frame(canvas)
        self.detailsFrame.grid(row=0, column=0)
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
    