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
from protlib_gui_ext import ToolTip, MultiListbox, centerWindows, askYesNo, configDefaults,\
    showInfo
from config_protocols import protDict, sections
from config_protocols import FontName, FontSize
from protlib_base import getProtocolFromModule, XmippProject,\
    getWorkingDirFromRunName, getExtRunName
from protlib_utils import reportError, runImageJPlugin, getProcessFromPid,\
    getRelatedProcess
from protlib_sql import SqliteDb, ProgramDb

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
ButtonSelectedColor = "DeepSkyBlue2"

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
        self.label = ProjectLabel(self, text=label_text)
        self.label.grid(row=0, column=0, sticky='sw')
        self.frameButtons = tk.Frame(self)
        self.frameButtons.grid(row=0, column=1, sticky='e')
        self.frameContent = tk.Frame(self)
        self.frameContent.grid(row=1, column=0, columnspan=2, sticky='nsew')
        
    def addButton(self, text, imagePath=None, **opts):
        btn = ProjectButton(self.frameButtons, text, imagePath, **opts)
        btn.pack(side=tk.LEFT, padx=5)
        return btn
    
    def displayButtons(self, show=True):
        if show:
            self.frameButtons.grid()
        else:
            self.frameButtons.grid_remove()
            
class ProjectButtonMenu(tk.Frame):
    def __init__(self, master, label_text, **opts):
        tk.Frame.__init__(self, master, **opts)
        self.label = ProjectLabel(self, text=label_text)
        self.label.pack(padx=5, pady=(5, 0))
        
    def addButton(self, button_text, **opts):
        btn = ProjectButton(self, button_text, **opts)
        btn.pack(fill=tk.X, padx=5, pady=(5, 0))
        return btn
        
class XmippProjectGUI():  
    def __init__(self, project):
        self.project = project
        
    def cleanProject(self):
        if askYesNo("DELETE confirmation", "You are going to DELETE all project data (runs, logs, results...Do you want to continue?", self.root):
            self.project.clean()
            self.close()
            
    def deleteTmpFiles(self):
        try:
            self.project.deleteTmpFiles()
            tkMessageBox.showinfo("Operation success", "All temporary files have been successfully removed")
        except Exception, e:
            tkMessageBox.showerror("Operation error ", str(e))
    
    def browseFiles(self):
        runImageJPlugin("512m", "XmippBrowser.txt", "", True)
        
    def initVariables(self):
        self.ToolbarButtonsDict = {}
        self.runButtonsDict = {}
        self.lastSelected = None
        self.lastRunSelected = None
        self.Frames = {}
        self.historyRefresh = None
        self.historyRefreshRate = 4 # Refresh rate will be 1, 2, 4 or 8 seconds
              
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
        self.fileMenu.add_command(label="Exit", command=self.close)
        #File menu
        self.menubar.add_cascade(label="File", menu=self.fileMenu)
        #Project menu
        self.menuProject = tk.Menu(self.root, tearoff=0)
        self.menubar.add_cascade(label="Project", menu=self.menuProject)
        self.menuProject.add_command(label="Browse files", command=self.browseFiles)
        self.menuProject.add_command(label="Remove temporary files", command=self.deleteTmpFiles)        
        self.menuProject.add_command(label="Clean project", command=self.cleanProject)
        
    def selectRunUpDown(self, event):
        if event.keycode == 111: # Up arrow
            self.lbHist.selection_move_up()
        elif event.keycode == 116: # Down arrow
            self.lbHist.selection_move_down()
        
    def createToolbarMenu(self, parent, opts=[]):
        if len(opts) > 0:
            menu = tk.Menu(parent, bg=ButtonBgColor, activebackground=ButtonBgColor, font=Fonts['button'], tearoff=0)
            prots = [protDict[o] for o in opts]
            for p in prots:
                #Following is a bit tricky, its due Python notion of scope, a for does not define a new scope
                # and menu items command setting
                def item_command(prot): 
                    def new_command(): 
                        self.launchProtocolGUI(self.project.newProtocol(prot.name))
                    return new_command 
                menu.add_command(label=p.title, command=item_command(p))
            menu.bind("<Leave>", self.unpostMenu)
            return menu
        return None

    
    #GUI for launching Xmipp Programs as Protocols        
    def launchProgramsGUI(self, event=None):
        text = protDict.xmipp_program.title
        last = self.lastSelected
        self.selectToolbarButton(text, False)
        if text != last and len(self.runs)>0:
            return
        db = ProgramDb()        
        root = tk.Toplevel()
        root.withdraw()
        root.title(text)
        root.columnconfigure(1, weight=1)
        root.rowconfigure(0, weight=1)
        root.rowconfigure(1, weight=1)
        detailsSection = ProjectSection(root, 'Details')
        txt = tk.Text(detailsSection.frameContent, width=50, height=10,
                        bg=BgColor, bd=1, relief=tk.RIDGE, font=Fonts['normal'])
        txt.pack(fill=tk.BOTH)
        detailsSection.grid(row=1, column=1, sticky='nsew', 
                            padx=5, pady=5)
        #Create programs panel
        progSection = ProjectSection(root, 'Programs')
        lb = tk.Listbox(progSection.frameContent, width=50, height=14,
                        bg=BgColor, bd=1, relief=tk.RIDGE, font=Fonts['normal'])

        def runClick(event=None):
            program_name = lb.get(int(lb.curselection()[0]))
            tmp_script = 'protocol_program_header.py'
            os.system(program_name + " --xmipp_write_protocol %(tmp_script)s " % locals())
            self.launchProtocolGUI(self.project.createRunFromScript(protDict.xmipp_program.name, 
                                                                    tmp_script, program_name))            
        
        def showSelection(event):
            #lb = event.widget
            program_name = lb.get(int(lb.curselection()[0]))
            program = db.selectProgram(program_name)
            txt.delete(1.0, tk.END)
            desc = program['usage'].splitlines()
            for line in desc:
                txt.insert(tk.END, line)
                
        detailsSection.addButton('Run', command=runClick)
        lb.bind('<ButtonRelease-1>', showSelection)
        lb.bind('<Double-Button-1>', runClick)
        lb.pack(fill=tk.BOTH)
        progSection.grid(row=0, column=1, sticky='nsew', padx=5, pady=5)
        #Create toolbar
        leftFrame = tk.Frame(root)
        leftFrame.grid(row=0, column=0, rowspan=2, sticky='nw', padx=5, pady=5)
        toolbar = tk.Frame(leftFrame, bd=2, relief=tk.RIDGE)
        section = ProjectButtonMenu(toolbar, 'Categories')
        section.pack(fill=tk.X, pady=5)
        categories = db.selectCategories()
        
        def fillListBox(programs):
            lb.delete(0, tk.END)
            for p in programs:
                lb.insert(tk.END, p['name'])
            
        def searchByKeywords(event=None):
            from protlib_xmipp import ProgramKeywordsRank
            keywords = searchVar.get().split()
            progRank = ProgramKeywordsRank(keywords)
            programs = db.selectPrograms()
            # Calculate ranks
            results = []
            for p in programs:
                rank = progRank.getRank(p)
                #Order by insertion sort
                if rank > 0:
                    name = p['name']
                    #self.maxlen = max(self.maxlen, len(name))
                    pos = len(results)
                    for i, e in reversed(list(enumerate(results))):
                        if e['rank'] < rank:
                            pos = i
                        else: break
                    results.insert(pos, {'rank': rank, 'name': name, 'usage': p['usage']})
            fillListBox(results)
            
        for c in categories:
            def addButton(c):
                section.addButton(c['name'], command=lambda:fillListBox(db.selectPrograms(c)))
            addButton(c)  
        toolbar.grid(row=0)
        ProjectLabel(leftFrame, text="Search").grid(row=1, pady=(10,5))
        searchVar = tk.StringVar()
        searchEntry = tk.Entry(leftFrame, textvariable=searchVar, bg=LabelBgColor)
        searchEntry.grid(row=2, sticky='ew')
        searchEntry.bind('<Return>', searchByKeywords)
        centerWindows(root, refWindows=self.root)
        root.deiconify()
        root.mainloop() 

    def launchRunJobMonitorGUI(self, run):
        runName = getExtRunName(run)
        text = run['script']
        root = tk.Toplevel()
        root.withdraw()
        root.title(text)
        root.columnconfigure(1, weight=1)
        root.rowconfigure(0, weight=1)
        root.rowconfigure(1, weight=1)
        
        def stopRun():
            if askYesNo("Confirm action", "Are you sure to stop run execution?" , parent=root):
                p = getProcessFromPid(run['pid'])
                p.terminateTree()
                self.project.projectDb.updateRunState(SqliteDb.RUN_ABORTED, run['run_id'])
                root.destroy()
                self.historyRefreshRate = 1
                self.updateRunHistory(self.lastSelected)
        
        detailsSection = ProjectSection(root, 'Process monitor')
        detailsSection.addButton("Stop run", command=stopRun)
        txt = tk.Text(detailsSection.frameContent, width=80, height=15,
                        bg=BgColor, bd=1, relief=tk.RIDGE, font=Fonts['normal'])
        txt.pack(fill=tk.BOTH)
        detailsSection.grid(row=1, column=1, sticky='nsew', padx=5, pady=5)
        
        def refreshInfo():
            p = getProcessFromPid(run['pid'])
            line = "Process id    : %(pid)s\nElapsed time  : %(etime)s\n\nSubprocess:\n" % p.info
            txt.delete(1.0, tk.END)
            txt.insert(tk.END, line)
            txt.insert(tk.END, "PID\t ARGS\t CPU(%)\t MEM(%)\n")
            childs = getRelatedProcess(getWorkingDirFromRunName(runName))
            for c in childs:
                c.info['pname'] = p.args.split()[0]
                line = "%(pid)s\t %(args)s\t %(pcpu)s\t %(pmem)s\n" % c.info
                txt.insert(tk.END, line)
            txt.after(3000, refreshInfo)
            
        refreshInfo()
        centerWindows(root, refWindows=self.root)
        root.deiconify()
        root.mainloop()
        
        
    def launchProtocolGUI(self, run, visualizeMode=False):
        run['group_name'] = self.lastSelected
        top = tk.Toplevel()
        gui = ProtocolGUI()
        gui.createGUI(self.project, run, top, 
                      lambda: self.protocolSaveCallback(run), visualizeMode)
        gui.launchGUI()
        
    def protocolSaveCallback(self, run):
        self.selectToolbarButton(run['group_name'], False)
        if self.lastSelected == run['group_name']:
            self.historyRefreshRate = 1
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
        if protGroup == 'All':
            self.runs = self.project.projectDb.selectRuns()
        else:
            self.runs = self.project.projectDb.selectRuns(protGroup)
        if len(self.runs) > 0:
            for run in self.runs:
                state = run['run_state']
                stateStr = SqliteDb.StateNames[state]
                if state == SqliteDb.RUN_STARTED:
                    stateStr += " - %d/%d" % self.project.projectDb.getRunProgress(run)
                self.lbHist.insert(tk.END, (getExtRunName(run), stateStr, run['last_modified']))   
            self.lbHist.selection_set(index)
            #Generate an automatic refresh after x ms, with x been 1, 2, 4 or 8 ms 
            self.historyRefresh = self.lbHist.after(self.historyRefreshRate*1000, self.updateRunHistory, protGroup, False)
            self.historyRefreshRate = min(2*self.historyRefreshRate, 8)
        else:
            self.updateRunSelection(-1)

    #---------------- Functions related with Popup menu ----------------------   
    def lastPair(self):
        if self.lastSelected:
            return  self.ToolbarButtonsDict[self.lastSelected]
        return None
        
    def unpostMenu(self, event=None):
        try: #I'm getting here an weird Tk exception
            menu = self.lastPair()[1]
            menu.unpost()
        except Exception:
            pass
            
    def postMenu(self, btn, menu):
        x, y, w = btn.winfo_x(), btn.winfo_y(), btn.winfo_width()
        xroot, yroot = self.root.winfo_x() + btn.master.winfo_x(), self.root.winfo_y()+ btn.master.winfo_y()
        menu.post(xroot + x + w + 10, yroot + y)
        
    def selectToolbarButton(self, text, showMenu=True):
        btn, menu = self.ToolbarButtonsDict[text]

        if self.lastSelected and self.lastSelected != text:
            lastBtn, lastMenu = self.lastPair()
            lastBtn.config(bg=ButtonBgColor, activebackground=ButtonActiveBgColor)
            if lastMenu:
                lastMenu.unpost()            
        
        if self.lastSelected != text:
            self.project.config.set('project', 'lastselected', text)
            self.project.writeConfig()
            self.updateRunHistory(text)            
            self.lastSelected = text  
            btn.config(bg=ButtonSelectedColor, activebackground=ButtonSelectedColor)
            
        if self.lastSelected and showMenu:
            self.postMenu(btn, menu)
            
    def updateRunSelection(self, index):
        state = tk.NORMAL
        details = self.Frames['details']
        if index == -1:
            state = tk.DISABLED
            #Hide details
            details.grid_remove()
        else:
            state = tk.NORMAL
            run = self.lastRunSelected
            prot = getProtocolFromModule(run['script'], self.project)
            if os.path.exists(prot.WorkingDir):
                summary = '\n'.join(prot.summary())
                showButtons = True
            else:
                summary = "This protocol run has not been executed yet"
                showButtons = False
            labels = [('Run', getExtRunName(run)),
                      ('\nCreated', run['init']),('   Modified', run['last_modified']), 
                      ('\nScript ', run['script']),('\nDirectory', prot.WorkingDir),
                      ('\nSummary', "\n"+summary)      ]
            self.detailsText.config(state=tk.NORMAL)
            self.detailsText.delete(1.0, tk.END)
            for k, v in labels:
                self.detailsText.insert(tk.END, '%s:   ' % k, 'tag_bold')
                self.detailsText.insert(tk.END, '%s' % v, 'tag_normal')
            self.detailsText.config(state=tk.DISABLED)
            details.displayButtons(showButtons)
            details.grid()
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
        state = run['run_state']
        if event == 'Edit':
            if state == SqliteDb.RUN_STARTED:
                self.launchRunJobMonitorGUI(run)
            else:
                self.launchProtocolGUI(run)
            #self.launchRunJobMonitorGUI(run)
        elif event == 'Copy':
            self.launchProtocolGUI(self.project.copyProtocol(run['protocol_name'], run['script']))
        elif event == "Delete":
            if askYesNo("Confirm DELETE", "All data related to this run will be DELETED. Do you want to continue?", self.root):
                self.project.deleteRun(run)
                self.updateRunHistory(self.lastSelected)
        elif event == "Visualize":
            pass
        
    def createToolbarFrame(self, parent):
        #Configure toolbar frame
        toolbar = tk.Frame(parent, bd=2, relief=tk.RIDGE)
        self.Frames['toolbar'] = toolbar
        #Create toolbar buttons
        #i = 1
        section = None
        for k, v in sections:
            section = ProjectButtonMenu(toolbar, k)
            section.pack(fill=tk.X, pady=5)
            for o in v:
                text = o[0]
                opts = o[1:]
                def btn_command(text): 
                    def new_command(): 
                        self.selectToolbarButton(text)
                    return new_command 
                btn = section.addButton(text, command=btn_command(text))
                menu = self.createToolbarMenu(section, opts)
                self.ToolbarButtonsDict[text] = (btn, menu)
        text = protDict.xmipp_program.title
        btn = section.addButton(text, command=self.launchProgramsGUI)
        self.ToolbarButtonsDict[text] = (btn, None)
        text = "All"
        btn = section.addButton(text, command=lambda: self.selectToolbarButton(text, False))
        self.ToolbarButtonsDict[text] = (btn, None)        
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
        def setupButton(k, v):
            btn =  history.addButton(k, v, command=lambda:self.runButtonClick(k))
            ToolTip(btn, k, 500)
            self.runButtonsDict[k] = btn
        for k, v in list:
            setupButton(k, v)
        self.lbHist = MultiListbox(history.frameContent, (('Run', 35), ('State', 15), ('Modified', 15)))
        self.lbHist.SelectCallback = self.runSelectCallback
        self.lbHist.DoubleClickCallback = lambda:self.runButtonClick("Edit")
        self.lbHist.AllowSort = False   
        return history     
        
    def createDetailsFrame(self, parent):
        details = ProjectSection(parent, 'Details')
        self.Frames['details'] = details
        #Create RUN details
        details.addButton("Analyze results", command=self.visualizeRun)
        details.addButton("Show output", command=self.showOutput)
        content = details.frameContent
        content.config(bg=BgColor, bd=1, relief=tk.RIDGE)
        content.grid_configure(pady=(5, 0))
        self.detailsText = tk.Text(content, height=15, width=70, border=0, 
                                   background='white', fg="black")
        self.detailsText.pack(fill=tk.BOTH)
        self.detailsText.tag_config('tag_normal', font=tkFont.Font(family=FontName, size=FontSize))
        self.detailsText.tag_config('tag_bold', font=tkFont.Font(family=FontName, size=FontSize, weight=tkFont.BOLD))
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
            self.selectToolbarButton(self.project.config.get('project', 'lastselected'), False)
        else:
            self.Frames['details'].grid_remove()
    
    def launchGUI(self, center=True):
        if center:
            centerWindows(self.root)
        self.root.deiconify()
        self.root.mainloop()
       
    def close(self, event=""):
        self.root.destroy()

    def showOutput(self, event=''):
        prot = getProtocolFromModule(self.lastRunSelected['script'], self.project)
        root = tk.Toplevel()
        root.withdraw()
        root.title("Output Console - %s" % self.lastRunSelected['script'])
        from protlib_gui_ext import OutputTextArea
        l = OutputTextArea(root, prot.LogPrefix)
        l.pack(side=tk.TOP)
        centerWindows(root, refWindows=self.root)
        root.deiconify()
        root.mainloop() 
        
    def visualizeRun(self, event=''):
        run = self.getLastRunDict()
        self.launchProtocolGUI(run, True)

from protlib_xmipp import XmippScript

class ScriptProtocols(XmippScript):
    def __init__(self):
        XmippScript.__init__(self, True)
        
    def defineParams(self):
        self.addUsageLine("Create Xmipp project on this folder.");
        ## params
        self.addParamsLine("[ -c  ]            : Clean project");
        self.addParamsLine("   alias --clean;");    
    
    def run(self):
        dir = os.getcwd()
        project = XmippProject(dir)
        launch = True
        if self.checkParam('--clean'):
            project.clean()
        else: #lauch project     
            if not project.exists():
                print 'You are in directory: ', dir
                answer = raw_input('Do you want to create a new xmipp_protocols PROJECT in this folder? [Y/n]:')
                if not answer or answer.lower() == 'y':
                    project.create()
                else: 
                    launch = False
            else:
                project.load()
        if launch:
            gui = XmippProjectGUI(project)
            gui.createGUI()
            gui.launchGUI()

if __name__ == '__main__':
    ScriptProtocols().tryRun()
