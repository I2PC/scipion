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
from protlib_gui import ProtocolGUI, Fonts, registerFont, registerCommonFonts
from protlib_gui_ext import ToolTip, MultiListbox, centerWindows
from config_protocols import protDict, sections
from protlib_base import getProtocolFromModule, XmippProject
from protlib_utils import reportError
from protlib_sql import SqliteDb

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

def configDefaults(opts, defaults):
    for key in defaults.keys():
        if not opts.has_key(key):
            opts[key] = defaults[key]
            
def ProjectButton(master, text, imagePath=None, **opts):
    configDefaults(opts, {'activebackground': ButtonActiveBgColor})
    btnImage = None
    if imagePath:
        try:
            from protlib_filesystem import getXmippPath
            imgPath = os.path.join(getXmippPath('resources'), imagePath)
            btnImage = tk.PhotoImage(file=imgPath)
        except tk.TclError:
            pass
    
    if btnImage:
        btn = tk.Button(master, image=btnImage, bd=0, height=28, width=28, **opts)
        btn.image = btnImage
    else:
        btn = tk.Button(master, text=text, font=Fonts['button'], bg=ButtonBgColor, **opts)
    return btn

    
def ProjectLabel(master, **opts):
    return tk.Label(master, font=Fonts['label'], fg=SectionTextColor, **opts)
        
class ProjectSection(tk.Frame):
    def __init__(self, master, label_text, **opts):
        tk.Frame.__init__(self, master, bd=2)
        self.config(**opts)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)
        self.label = tk.Label(self, text=label_text, font=Fonts['label'], fg=SectionTextColor)
        self.label.grid(row=0, column=0, sticky='sw')
        self.frameButtons = tk.Frame(self)
        self.frameButtons.grid(row=0, column=1, sticky='e')
        self.frameContent = tk.Frame(self)
        self.frameContent.grid(row=1, column=0, columnspan=2, sticky='nsew')
        
    def addButton(self, text, imagePath=None, **opts):
        btn = ProjectButton(self.frameButtons, text, imagePath, **opts)
        btn.pack(side=tk.LEFT, padx=5)
        return btn
        
class XmippProjectGUI():  
    def __init__(self, project):
        self.project = project
        
    def deleteProject(self):
        if tkMessageBox.askyesno("DELETE confirmation", "You are going to DELETE all project data (runs, logs, results...Do you want to continue?"):
            self.project.clean()
            self.close()
      
    def initVariables(self):
        self.ToolbarButtonsDict = {}
        self.runButtonsDict = {}
        self.lastSelected = None
        self.lastRunSelected = None
        self.Frames = {}
        self.historyRefresh = None
              
    def addBindings(self):
        self.root.bind('<Configure>', self.unpostMenu)
        self.root.bind("<Unmap>", self.unpostMenu)
        self.root.bind("<Map>", self.unpostMenu)
        self.root.bind('<Return>', lambda e: self.runButtonClick('Edit'))
        self.root.bind('<Control_L><Return>', lambda e: self.runButtonClick('Copy'))
        self.root.bind('<Delete>', lambda e: self.runButtonClick('Delete'))
        self.root.bind('<Up>', self.selectRunUpDown)
        self.root.bind('<Down>', self.selectRunUpDown)
        self.root.bind('<Alt_L><c>', self.close )
        self.root.bind('<Alt_L><l>', self.showOutput)
        self.root.bind('<Alt_L><a>', self.visualizeRun)
        
    def createMainMenu(self):
        self.menubar = tk.Menu(self.root)
        self.fileMenu = tk.Menu(self.root, tearoff=0)
        self.fileMenu.add_command(label="Delete project", command=self.deleteProject)
        self.fileMenu.add_command(label="Exit", command=self.close)
        self.menubar.add_cascade(label="File", menu=self.fileMenu)
        
    def selectRunUpDown(self, event):
        if event.keycode == 111: # Up arrow
            self.lbHist.selection_move_up()
        elif event.keycode == 116: # Down arrow
            self.lbHist.selection_move_down()
        
    def createToolbarButton(self, row, text, opts=[]):
        '''Add a button to left toolbar'''
        btn = tk.Button(self.Frames['toolbar'], bd = 1, text=text, font=Fonts['button'], relief=tk.RAISED,
                         bg=ButtonBgColor, activebackground=ButtonBgColor)
        btn.grid(row = row, column = 0, sticky='ew', pady=2, padx=5)
        if len(opts) > 0:
            menu = tk.Menu(self.root, bg=ButtonBgColor, activebackground=ButtonBgColor, font=Fonts['button'], tearoff=0)
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

    def launchProtocolGUI(self, run, visualizeMode=False):
        run['group_name'] = self.lastSelected
        top = tk.Toplevel()
        gui = ProtocolGUI()
        gui.createGUI(self.project, run, top, 
                      lambda: self.protocolSaveCallback(run), visualizeMode)
        gui.fillGUI()
        gui.launchGUI()
        
    def protocolSaveCallback(self, run):
        self.selectToolbarButton(run['group_name'], False)
        if self.lastSelected == run['group_name']:
            self.updateRunHistory(self.lastSelected)
          
        
    def updateRunHistory(self, protGroup, selectFirst=True):
        #Cancel if there are pending refresh
        if self.historyRefresh:
            self.lbHist.after_cancel(self.historyRefresh)
            self.historyRefresh = None
        index = 0
        if not selectFirst:
            index = self.lbHist.selectedIndex()
        self.lbHist.delete(0, tk.END)
        self.runs = self.project.projectDb.selectRuns(protGroup)
        if len(self.runs) > 0:
            for run in self.runs:
                run_name = '%s_%s' % (run['protocol_name'], run['run_name'])
                state = run['run_state']
                stateStr = SqliteDb.StateNames[state]
                if state == SqliteDb.RUN_STARTED:
                    stateStr += " - %d/%d" % self.project.projectDb.getRunProgress(run)
                self.lbHist.insert(tk.END, (run_name, stateStr, run['last_modified']))   
            self.lbHist.selection_set(index)
            #Generate an automatic refresh after 3000 ms 
            self.historyRefresh = self.lbHist.after(1000, self.updateRunHistory, protGroup, False)
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
        
    def selectToolbarButton(self, text, showMenu=True):
        btn, menu = self.ToolbarButtonsDict[text]

        if self.lastSelected and self.lastSelected != text:
            lastBtn, lastMenu = self.lastPair()
            lastBtn.config(bg=ButtonBgColor)
            lastMenu.unpost()            
        
        if self.lastSelected and showMenu:
            self.postMenu(btn, menu)
        
        if self.lastSelected != text:
            self.project.config.set('project', 'lastselected', text)
            self.project.writeConfig()
            self.updateRunHistory(text)            
        self.lastSelected = text  
            
    def updateRunSelection(self, index):
        state = tk.NORMAL
        if index == -1:
            state = tk.DISABLED
            #Hide details
            self.Frames['details'].grid_remove()
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
            self.Frames['details'].grid()
        for btn in self.runButtonsDict.values():
            btn.config(state=state)
            
    def runSelectCallback(self, index):
        if index >= 0:
            self.lastRunSelected = self.runs[index]
        self.updateRunSelection(index)
            
    def getLastRunDict(self):
        from protlib_sql import runColumns
        run = dict(zip(runColumns, self.lastRunSelected))
        run['source'] = run['script']        
        return run
    
    def runButtonClick(self, event=None):
        run = self.getLastRunDict()
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
        

    def createToolbarFrame(self, parent):
        #Configure toolbar frame
        toolbar = tk.Frame(parent, bd=2, relief=tk.RIDGE)
        self.Frames['toolbar'] = toolbar
        #Create toolbar buttons
        i = 1
        for k, v in sections:
            ProjectLabel(toolbar, text=k).grid(column=0, row=i)
            i += 1
            for btn in v:
                self.createToolbarButton(i, btn[0], btn[1:])
                i += 1
        return toolbar
                
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
            btn = tk.Button(frame, text=text, font=Fonts['button'], bg=ButtonBgColor)
        btn.config(command=lambda:self.runButtonClick(text), 
                 activebackground=ButtonActiveBgColor)
        btn.grid(row=0, column=col)
        ToolTip(btn, text, 500)
        self.runButtonsDict[text] = btn
    
    def createHistoryFrame(self, parent):
        history = ProjectSection(parent, 'History')
        self.Frames['history'] = history
        list = [('Edit', 'edit.gif'), ('Copy', 'copy.gif'), ('Delete', 'delete.gif')]
        for k, v in list:
            btn =  history.addButton(k, v, command=lambda:self.runButtonClick(k))
            ToolTip(btn, k, 500)
            self.runButtonsDict[k] = btn
            
        self.lbHist = MultiListbox(history.frameContent, (('Run', 35), ('State', 15), ('Modified', 15)))
        self.lbHist.SelectCallback = self.runSelectCallback
        self.lbHist.DoubleClickCallback = lambda:self.runButtonClick("Edit")
        self.lbHist.AllowSort = False   
        return history     
        
    def addDetailsLabel(self, parent, text, row, col, colspan=1):
        label = tk.Label(parent, text=text, font=Fonts['details'], bg=BgColor)
        label.grid(row=row, column=col, sticky='ne', padx=5)
        label = tk.Label(parent,text="", font=Fonts['details_bold'],
                      bg=BgColor, justify=tk.LEFT)
        label.grid(row=row, column=col+1, sticky='nw', padx=5, columnspan=colspan)
        self.DetailsLabelsDict[text] = label
        
    def createDetailsFrame(self, parent):
        details = ProjectSection(parent, 'Details')
        self.Frames['details'] = details
        #Create RUN details
        details.addButton("Analyze results", command=self.visualizeRun)
        details.addButton("Show output", command=self.showOutput)
        content = details.frameContent
        content.config(bg=BgColor, bd=1, relief=tk.RIDGE)
        content.grid_configure(pady=(5, 0))
        self.DetailsLabelsDict = {}
        registerFont('details', family=FontName, size=FontSize-1, weight=tkFont.BOLD)
        registerFont('details_bold', family=FontName, size=FontSize-1)
        self.addDetailsLabel(content, 'Run:', 0, 0)
        self.addDetailsLabel(content, 'Protocol:', 1, 0)
        self.addDetailsLabel(content, 'Created:', 0, 2)
        self.addDetailsLabel(content, 'Modified:', 1, 2)
        self.addDetailsLabel(content, 'Summary:', 2, 0, 3)
        return details

    def createGUI(self, root=None):
        if not root:
            root = tk.Tk()
        self.root = root
        root.withdraw() # Hide the windows for centering
        self.root.title("Xmipp Protocols")
        self.initVariables()        
        self.createMainMenu()
        self.addBindings()
        #Create main frame that will contain all other sections
        #Configure min size and expanding behaviour
        root.minsize(750, 500)
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        main = tk.Frame(self.root)
        main.grid(row=0, column=0, sticky="nsew")
        main.columnconfigure(0, minsize=150)
        main.columnconfigure(1, minsize=450, weight=1)
        main.rowconfigure(1, minsize=200, weight=1)
        main.rowconfigure(2, minsize=200, weight=1)
        self.Frames['main'] = main
        
        registerCommonFonts()
        
        #Create section frames and locate them
        self.createToolbarFrame(main).grid(row=1, column=0, sticky='nse', padx=5, pady=5, rowspan=2)
        self.createHistoryFrame(main).grid(row=1, column=1, sticky='nsew', padx=5, pady=5)
        self.createDetailsFrame(main).grid(row=2, column=1, sticky='nsew', padx=5, pady=5)
        
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
       
    def close(self, event=""):
        self.root.destroy()

    def showOutput(self, event=''):
        prot = getProtocolFromModule(self.lastRunSelected['script'], self.project)
        root = tk.Toplevel()
        root.title("Output Console - %s" % self.lastRunSelected['script'])
        from protlib_gui_ext import OutputTextArea
        l = OutputTextArea(root, prot.LogPrefix)
        l.pack(side=tk.TOP)
        root.mainloop() 
        
    def visualizeRun(self, event=''):
        run = self.getLastRunDict()
        self.launchProtocolGUI(run, True)

if __name__ == '__main__':
    dir = os.getcwd()
    project = XmippProject(dir)    
    import sys
    if len(sys.argv) > 1:
        if sys.argv[1] == "--clean":
            project.clean()
        else:
            reportError('Unrecognized option: %s' % sys.argv[1])
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
    
