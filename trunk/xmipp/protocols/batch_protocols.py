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
import Tkinter as tk
import tkMessageBox
import tkFont
from protlib_gui import ProtocolGUI
from protlib_gui_ext import ToolTip, MultiListbox, centerWindows
from config_protocols import protDict, sections
from protlib_base import getProtocolFromModule, XmippProject

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


        
class XmippProjectGUI():  
    def __init__(self, project):
        self.project = project
        
    def deleteProject(self):
        if tkMessageBox.askyesno("DELETE confirmation", "You are going to DELETE all project data (runs, logs, results...Do you want to continue?"):
            self.project.clean()
            self.onExit()            
            
    def createMainMenu(self):
        self.menubar = tk.Menu(self.root)
        self.fileMenu = tk.Menu(self.root, tearoff=0)
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
        self.root.bind('<Return>', lambda e: self.runButtonClick('Edit'))
        self.root.bind('<Delete>', lambda e: self.runButtonClick('Delete'))
        self.root.bind('<Up>', self.selectRunUpDown)
        self.root.bind('<Down>', self.selectRunUpDown)
        
    def selectRunUpDown(self, event):
        if event.keycode == 111: # Up arrow
            self.lbHist.selection_move_up()
        elif event.keycode == 116: # Down arrow
            self.lbHist.selection_move_down()
            
    def addHeaderLabel(self, parent, text, row, col=0):
        '''Add a label to left toolbar'''
        label = tk.Label(parent, text=text, font=self.LabelFont, fg=SectionTextColor)
        label.grid(row = row, column=col)
        return label
        
    def createToolbarButton(self, row, text, opts=[]):
        '''Add a button to left toolbar'''
        btn = tk.Button(self.toolbar, bd = 1, text=text, font=self.ButtonFont, relief=tk.RAISED,
                         bg=ButtonBgColor, activebackground=ButtonBgColor)
        btn.grid(row = row, column = 0, sticky='ew', pady=2, padx=5)
        if len(opts) > 0:
            menu = tk.Menu(self.root, bg=ButtonBgColor, activebackground=ButtonBgColor, font=self.ButtonFont, tearoff=0)
            prots = [protDict.protocolDict[o] for o in opts]
            for p in prots:
                #Following is a bit tricky, its due Python notion of scope, a for does not define a new scope
                # and menu items command setting
                def item_command(prot): 
                    def new_command(): 
                        self.launchProtocolGUI(self.project.newProtocol(prot.key))
                    return new_command 
                menu.add_command(label=p.title, command=item_command(p))
            menu.bind("<Leave>", self.unpostMenu)
        self.ToolbarButtonsDict[text] = (btn, menu)
        btn.config(command=lambda:self.selectToolbarButton(text))

    def launchProtocolGUI(self, run):
        run['group_name'] = self.lastSelected
        top = tk.Toplevel()
        gui = ProtocolGUI()
        gui.createGUI(self.project, run, top, lambda: self.protocolSaveCallback(run))
        gui.fillGUI()
        gui.launchGUI()
        
    def protocolSaveCallback(self, run):
        if self.lastSelected == run['group_name']:
            self.updateRunHistory(self.lastSelected) 

    def updateRunHistory(self, protGroup): 
        self.runs = self.project.projectDb.selectRuns(protGroup)
        self.lbHist.delete(0, tk.END)
        for run in self.runs:
            self.lbHist.insert(tk.END, ('%s_%s' % (run['protocol_name'], run['run_name']),
                                run['last_modified']))   
        if len(self.runs) > 0:
            self.lbHist.selection_set(0)
        else:
            self.updateRunSelection(-1)

    #---------------- Functions related with Popup menu ----------------------   
    def lastPair(self):
        if self.lastSelected:
            return  self.ToolbarButtonsDict[self.lastSelected]
        return None
        
    def unpostMenu(self, event=None):
        if self.lastSelected:
            menu = self.lastPair()[1]
            if menu:
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
            print 'posting menu'
            self.postMenu(btn, menu)
        
        if self.lastSelected != text:
            self.project.config.set('project', 'lastselected', text)
            self.project.writeConfig()
            #self.updateRunSelection(-1)
            print 'updating'
            self.updateRunHistory(text)            
        self.lastSelected = text  
            
    def updateRunSelection(self, index):
        state = tk.NORMAL
        if index == -1:
            state = tk.DISABLED
            #Hide details
            self.frameDetails.grid_remove()
            self.buttonDetails.grid_remove()
        else:
            state = tk.NORMAL
            #Show details
            run = self.lastRunSelected
            self.DetailsLabelsDict['Run:'].config(text=run['run_name'])
            self.DetailsLabelsDict['Protocol:'].config(text=run['protocol_name'])
            self.DetailsLabelsDict['Created:'].config(text=run['init'])
            self.DetailsLabelsDict['Modified:'].config(text=run['last_modified'])
            prot = getProtocolFromModule(run['script'], self.project)
            summary = '\n'.join(prot.summary())
            self.DetailsLabelsDict['Summary:'].config(text=summary)
            self.frameDetails.grid(row=4, column=1,sticky='nsew', columnspan=2)
            self.buttonDetails.grid()
        for btn in self.runButtonsDict.values():
            btn.config(state=state)
            
    def runSelectCallback(self, index):
        if index >= 0:
            self.lastRunSelected = self.runs[index]
        self.updateRunSelection(index)
            
    def runButtonClick(self, event=None):
        from protlib_sql import runColumns
        run = dict(zip(runColumns, self.lastRunSelected))
        run['source'] = run['script']        
        if event == 'Edit':
            self.launchProtocolGUI(run)
        elif event == 'Copy':
            self.launchProtocolGUI(self.project.copyProtocol(run['protocol_name'], run['script']))
        elif event == "Delete":
            if tkMessageBox.askyesno("Confirm DELETE", "All data related to this run will be DELETED. Do you want to continue?"):
                self.project.deleteRun(run)
                self.updateRunHistory(self.lastSelected)
        elif event == "Visualize":
            pass
        

    def createToolbar(self):
        #Configure toolbar frame
        self.toolbar = tk.Frame(self.frame, bd=2, relief=tk.RIDGE)
        self.toolbar.grid(row=0, column=0, sticky='nws', 
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
                from protlib_filesystem import getXmippPath
                imgPath = os.path.join(getXmippPath('resources'), imageFilename)
                btnImage = tk.PhotoImage(file=imgPath)
            except tk.TclError:
                pass
        
        if btnImage:
            btn = tk.Button(frame, image=btnImage, bd=0, height=28, width=28)
            btn.image = btnImage
        else:
            btn = tk.Button(frame, text=text, font=self.ButtonFont, bg=ButtonBgColor)
        btn.config(command=lambda:self.runButtonClick(text), 
                 activebackground=ButtonActiveBgColor)
        btn.grid(row=0, column=col)
        ToolTip(btn, text, 500)
        self.runButtonsDict[text] = btn
    
    def createRunHistory(self):
        self.addHeaderLabel(self.frame, 'History', 0, 1)
        #Button(self.frame, text="Edit").grid(row=0, column=2)
        frame = tk.Frame(self.frame)
        frame.grid(row=0, column=2)
        self.addRunButton(frame, "Edit", 0, 'edit.gif')
        self.addRunButton(frame, "Copy", 1, 'copy.gif')
        #self.addRunButton(frame, "Visualize", 2, 'visualize.gif')
        self.addRunButton(frame, "Delete", 2, 'delete.gif')
        #self.addRunButton(frame, "Help", 4, 'help.gif')
        self.frameHist = tk.Frame(self.frame)
        self.frameHist.grid(row=2, column=1, sticky='nsew', columnspan=2, padx=5, pady=(0, 5))
        self.lbHist = MultiListbox(self.frameHist, (('Run', 40), ('Modified', 20)))
        self.lbHist.SelectCallback = self.runSelectCallback
        self.lbHist.DoubleClickCallback = lambda:self.runButtonClick("Edit")
        self.lbHist.AllowSort = False        
        #self.lbHist.pack()
        
    def addDetailsLabel(self, text, row, col, sumCol=True):
        label = tk.Label(self.frameDetails,text=text, font=self.DetailsFontBold, bg=BgColor)
        label.grid(row=row, column=col, sticky='ne', padx=5)
        if sumCol:
            col += 1
        else:
            row += 1
        label = tk.Label(self.frameDetails,text="", font=self.DetailsFont, 
                      bg=BgColor, justify=tk.LEFT)
        colspan = 1
        if text == 'Summary:':
            colspan = 3
        label.grid(row=row, column=col, sticky='nw', padx=5, columnspan=colspan)
        self.DetailsLabelsDict[text] = label
        
    def createRunDetails(self):
        #Prepare fonts
        self.DetailsFontBold = tkFont.Font(family=FontName, size=FontSize-1, weight=tkFont.BOLD)
        self.DetailsFont = tkFont.Font(family=FontName, size=FontSize-1)
        #Create RUN details
        self.addHeaderLabel(self.frame, 'Details', 3, 1)
        self.buttonDetails = tk.Button(self.frame, text="Analyze results", font=self.ButtonFont, 
                                    bg=ButtonBgColor, activebackground=ButtonActiveBgColor,
                                    command=self.visualizeRun)
        self.buttonDetails.grid(row=3, column=2, padx=5, pady=5)
        self.frameDetails = tk.Frame(self.frame, bg=BgColor, bd=1, relief=tk.RIDGE)
        self.frameDetails.grid(row=4, column=1, sticky='nsew', columnspan=2, padx=5, pady=5)
        self.DetailsLabelsDict = {}
        self.addDetailsLabel('Run:', 0, 0)
        self.addDetailsLabel('Protocol:', 1, 0)
        self.addDetailsLabel('Created:', 0, 2)
        self.addDetailsLabel('Modified:', 1, 2)
        self.addDetailsLabel('Summary:', 2, 0)

    def createGUI(self, root=None):
        if not root:
            root = tk.Tk()
        self.root = root
        root.withdraw() # Hide the windows for centering
        self.root.title("Xmipp Protocols")
        self.createMainMenu()
        #Create a main frame that contains all other widgets
        self.frame = tk.Frame(self.root)
        self.frame.pack(fill=tk.BOTH)
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

    def visualizeRun(self):
        getProtocolFromModule(self.lastRunSelected['script'], self.project).visualize()

if __name__ == '__main__':
    dir = os.getcwd()
    project = XmippProject(dir)    
    import sys
    if len(sys.argv) > 1:
        # Launch a protocol directly
        script = sys.argv[1]
        project.load()  
        gui = ProtocolGUI()
        gui.createGUI(script)
        gui.launchGUI()
     
    else: #lauch project     
        if not project.exists():
            print 'You are in directory: ', dir
            answer = raw_input('Do you want to create a new xmipp_protocols PROJECT in this folder? [Y/n]:')
            if not answer or answer.lower() == 'y':
                project.create()
        else:
            project.load()
        gui = XmippProjectGUI(project)
        gui.createGUI()
        gui.launchGUI()
    
