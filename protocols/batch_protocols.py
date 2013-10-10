#!/usr/bin/env xmipp_python
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
import time
import threading
from os.path import exists
import Tkinter as tk
#import protlib_gui_figure
from protlib_gui import ProtocolGUI, Fonts, registerCommonFonts
from protlib_gui_ext import ToolTip, centerWindows, askYesNo, showInfo, XmippTree, \
    showBrowseDialog, showTextfileViewer, showError, TaggedText, XmippButton, ProjectLabel,\
    FlashMessage, showWarning, YesNoDialog, getXmippImage
from config_protocols import *
from protlib_base import XmippProject, getExtendedRunName, splitExtendedRunName,\
    getScriptFromRunName, ProtocolExecutor
from protlib_utils import ProcessManager,  getHostname, loadModule, datetime_fromdb, pretty_delta
from protlib_sql import SqliteDb, ProgramDb
from protlib_filesystem import getXmippPath, copyFile
from protlib_parser import ProtocolParser

# Redefine BgColor
BgColor = "white"

ACTION_EDIT = 'Edit'
ACTION_COPY = 'Copy'
ACTION_DELETE = 'Delete'
ACTION_REFRESH = 'Refresh'
ACTION_STEPS = 'Show_Steps'
ACTION_TREE = 'Show_DepsTree'
ACTION_STOP = 'Stop'
ACTION_DEFAULT = 'Default'

ActionIcons = {
    ACTION_EDIT:  'edit.gif',
    ACTION_COPY:  'copy.gif',
    ACTION_DELETE:  'delete.gif',
    ACTION_REFRESH:  'refresh.gif',
    ACTION_STEPS:  'run_steps.gif',
    ACTION_TREE:  'tree2.gif',
    ACTION_STOP: 'stop.gif'
               }


GROUP_ALL = 'All'
#GROUP_XMIPP = protDict.xmipp.title
        
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
        btn = XmippButton(self.frameButtons, text, imagePath, **opts)
        btn.pack(side=tk.LEFT, padx=5)
        return btn
    
    def addWidget(self, Widget, **opts):
        w = Widget(self.frameButtons, **opts)
        w.pack(side=tk.LEFT, padx=5)
        return w
    
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
        btn = XmippButton(self, button_text, **opts)
        btn.pack(fill=tk.X, padx=5, pady=(5, 0))
        return btn
        
class XmippProjectGUI():  
    def __init__(self, project):
        self.project = project
        #self.pm = ProcessManager()
        
    def cleanProject(self):
        if askYesNo("DELETE confirmation", "You are going to <DELETE ALL> project data (runs, logs, results...)\nDo you really want to continue?", self.root):
            self.project.clean()
            self.close()
            
    def deleteTmpFiles(self):
        try:
            self.project.deleteTmpFiles()
            showInfo("Operation success", "All temporary files have been successfully removed", self.root)
        except Exception, e:
            showError("Operation error ", str(e), self.root)
    
    def browseFiles(self):
        showBrowseDialog(parent=self.root, seltype="none", selmode="browse")
        
    def initVariables(self):
        self.ToolbarButtonsDict = {}
        self.runButtonsDict = {}
        self.lastDisplayGroup = None
        self.lastRunSelected = None
        self.lastMenu = None
        self.Frames = {}
        self.historyRefresh = None
        self.historyRefreshRate = 4 # Refresh rate will be 1, 2, 4 or 8 seconds
        ## Create some icons
        self.images = {}
        self.graphView = False
        ## This are images icons for each run_state
        paths = ['save_small', 'select_run', 'progress_none', 
                 'progress_ok', 'progress_error', 'level_warning']
        for i, path in enumerate(paths):
            self.images[i] = getXmippImage(path + '.gif')
            
    def getStateImage(self, state):
        return self.images[state]
    
    def getImage(self, imgName):
        ''' Create the image and store in dictionary'''
        if not imgName in self.images.keys():
            self.images[imgName] = getXmippImage(imgName + '.gif')
        return self.images[imgName]
              
    def addBindings(self):
        #self.root.bind('<FocusOut>', self.unpostMenu)
        self.root.bind("<Key>", self.unpostMenu)
        #self.root.bind('<Button-1>', self.unpostMenu)
        self.root.bind('<Return>', lambda e: self.runButtonClick(ACTION_DEFAULT))
        self.root.bind('<Control_L><Return>', lambda e: self.runButtonClick(ACTION_COPY))
        self.root.bind('<Alt_L><e>', lambda e: self.runButtonClick(ACTION_EDIT))
        self.root.bind('<Alt_L><d>', lambda e: self.runButtonClick(ACTION_COPY))
        self.root.bind('<Delete>', lambda e: self.runButtonClick(ACTION_DELETE))
        self.root.bind('<Alt_L><t>', lambda e: self.runButtonClick(ACTION_TREE))
        self.root.bind('<Alt_L><c>', self.close )
        self.root.bind('<Alt_L><o>', self.showOutput)
        self.root.bind('<Alt_L><a>', self.visualizeRun)
        self.root.bind('<Up>', self.treeHist.selection_up)
        self.root.bind('<Down>', self.treeHist.selection_down)
        self.root.bind("<MouseWheel>", self.graphHist.onScroll)
        
    def createMainMenu(self):
        self.menubar = tk.Menu(self.root)
        #Project menu
        self.menuProject = tk.Menu(self.root, tearoff=0)
        self.menubar.add_cascade(label="Project", menu=self.menuProject)
        self.browseFolderImg = tk.PhotoImage(file=getXmippPath('resources', 'folderopen.gif'))
        self.menuProject.add_command(label=" Browse files", command=self.browseFiles, 
                                     image=self.getImage('folderopen'), compound=tk.LEFT)
        self.menuProject.add_command(label=" Remove temporary files", command=self.deleteTmpFiles,
                                     image=self.getImage('delete'), compound=tk.LEFT) 
        self.menuProject.add_command(label=" Clean project", command=self.cleanProject,
                                     image=self.getImage('clean'), compound=tk.LEFT)
        self.menuProject.add_separator()
        self.menuProject.add_command(label=" Exit", command=self.close)
        # Help menu
        self.menuHelp = tk.Menu(self.root, tearoff=0)
        self.menubar.add_cascade(label="Help", menu=self.menuHelp)
        self.menuHelp.add_command(label=" Online help", command=self.openHelp,
                                     image=self.getImage('online_help'), compound=tk.LEFT)
        self.menuHelp.add_command(label=" About Xmipp", command=self.openAbout,
                                     compound=tk.LEFT)        
        self.menuRun = tk.Menu(self.root, tearoff=0)
        
    def openHelp(self, e=None):
        from protlib_gui_ext import openLink
        openLink('http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/WebHome')
        
    def openAbout(self, e=None):
        from about_gui import createAboutDialog
        createAboutDialog()
        
    def selectRunUpDown(self, event):
        if event.keycode == 111: # Up arrow
            self.treeHist.selection_move_up()
        elif event.keycode == 116: # Down arrow
            self.treeHist.selection_move_down()
        
    def createToolbarMenu(self, parent, opts=[]):
        bgColor = ButtonBgColor
        import sys
        if sys.platform  == 'darwin':
            # This is needed because in Mac the bgcolor is not working properly
            bgColor='black'
        if len(opts) > 0:
            menu = tk.Menu(parent, bg=bgColor, activebackground=bgColor, font=Fonts['button'], tearoff=0)
#            prots = [protDict[o] for o in opts]
            for o in opts:
                #Following is a bit tricky, its due Python notion of scope, a for does not define a new scope
                # and menu items command setting
                def item_command(prot): 
                    def new_command(): 
                        self.launchProtocolGUI(self.project.newProtocol(prot.name))
                    return new_command 
                # If p is a list, means is a submenu
                if type(o) == list:
                    menu2 = tk.Menu(menu, bg=bgColor, activebackground=bgColor, font=Fonts['button'], tearoff=0)
                    menu.add_cascade(label=o[0], menu=menu2)
                    for o2 in o[1:]:
                        p = protDict[o2]
                        menu2.add_command(label=p.title, command=item_command(p))
                    #menu2.bind("<Leave>", self.unpostMenu)
                else:
                    p = protDict[o]                  
                    # Treat special case of "custom" protocol
                    if p.title == 'Custom':
                        menu.add_command(label=p.title, command=self.selectCustomProtocol)
                    else:
                        menu.add_command(label=p.title, command=item_command(p))
            #menu.bind("<Leave>", self.unpostMenu)
            return menu
        return None

    def selectCustomProtocol(self, e=None):
        selection = showBrowseDialog(parent=self.root, seltype="file", filter="custom*py")
        if selection is not None:
            customProtFile = selection[0]
            mod = loadModule(customProtFile)
            prefix = getattr(mod, 'CustomPrefix', '')
            run = self.project.createRunFromScript(protDict.custom.name, customProtFile, prefix=prefix)
            #from os.path import basename
            #tmpFile = self.project.projectTmpPath(basename(run['script']))
            #copyFile(None, customProtFile, tmpFile)
            #run['source'] = tmpFile
            self.launchProtocolGUI(run)
        
    #GUI for launching Xmipp Programs as Protocols        
    def launchProgramsGUI(self, event=None):
        text = protDict.xmipp.title
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
            tmp_script = self.project.projectTmpPath('protocol_program_header.py')
            os.system(program_name + " --xmipp_write_protocol %(tmp_script)s " % locals())
            self.launchProtocolGUI(self.project.createRunFromScript(protDict.xmipp.name, 
                                                                    tmp_script, program_name))
            root.destroy();           
        
        def showSelection(event):
            #lb = event.widget
            program_name = lb.get(int(lb.curselection()[0]))
            program = db.selectProgram(program_name)
            txt.delete(1.0, tk.END)
            desc = program['usage'].splitlines()
            desc = desc[:5]# show max of 5 lines
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

    def stopRun(self):        
        if askYesNo("Confirm action", "Are you sure to <STOP> run execution?" , parent=self.root):
            run = self.getLastRunDict()
            pm = ProcessManager(run)
            pm.stopProcessGroup(self.project)
            
    def launchRunJobMonitorGUI(self, run):
        runName = getExtendedRunName(run)
        script = run['script']
        jobId = run['jobid']
        if jobId == SqliteDb.UNKNOWN_JOBID:
            showInfo('No jobId tracked', "The <jobId> of this process wasn't tracked", parent=self.root)
            return
        pm = ProcessManager(run)
        root = tk.Toplevel()
        root.withdraw()
        root.title(script)
        root.columnconfigure(1, weight=1)
        #root.rowconfigure(0, weight=1)
        root.rowconfigure(1, weight=1)
	
        def onClose():
            self.process_check = False 
            root.destroy()
	    
        def updateAndClose():
            root.destroy()
            self.historyRefreshRate = 1
            self.process_check = False
            self.updateRunHistory(self.lastDisplayGroup)
            
        def stopRun():
            if askYesNo("Confirm action", "Are you sure to <STOP> run execution?" , parent=root):
                pm.stopProcessGroup(self.project)
                updateAndClose()
                
        def processWorker():
            '''This function will be invoke from a worker thread 
            to update process info and avoid freeze the GUI '''
            self.process = pm.getProcessFromPid()
            while self.process_check and self.process:
                self.process_childs = pm.getProcessGroup()
                time.sleep(4)
                self.process = pm.getProcessFromPid()
                
        def refreshInfo():
            #p = pm.getProcessFromPid()
            #txt.delete(1.0, tk.END)
            tree.clear()
            if self.process:
                labPid.set(self.process.pid)
                labElap.set(self.process.etime)
                #line = "Process id    : %(pid)s\nElapsed time  : %(etime)s\n\nSubprocess:\n" % p.info
                #txt.insert(tk.END, line)
                #txt.insert(tk.END, "PID\t ARGS\t CPU(%)\t MEM(%)\n")
                #childs = pm.getProcessGroup()
                lastHost = ''
                for c in self.process_childs:
                    c.info['pname'] = pname = os.path.basename(c.args.split()[0])
                    if pname not in ['grep', 'python', 'sh', 'bash', 'xmipp_python', 'mpirun']:
                        #line = "%(pid)s\t %(pname)s\t %(pcpu)s\t %(pmem)s\n"
                        
                        if c.host != lastHost:
                            #txt.insert(tk.END, "%s\n" % c.host, 'bold')
                            tree.insert('', 'end', c.host, text=c.host)
                            tree.item(c.host, open=True)
                            lastHost = c.host
                        #txt.insert(tk.END, line % c.info, 'normal')
                        tree.insert(lastHost, 'end', values=(c.pid, c.pcpu, c.pmem, pname))
                tree.after(4500, refreshInfo)
            else:
                updateAndClose()
                
        detailsSection = ProjectSection(root, 'Process monitor')
        detailsSection.addButton("Stop run", command=stopRun)
        cols = ('pid', '%cpu', '%mem', 'command')
        frame = detailsSection.frameContent
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(2, weight=1)
        tree = XmippTree(frame, columns=cols)
        for c in cols:
            tree.heading(c, text=c)
            if c != 'command':
                tree.column(c, width=60)
        tree.heading('#0', text='Node')
        tk.Label(frame, text='Process id:', font=Fonts['button']).grid(row=0, column=0, sticky='e')
        tk.Label(frame, text='Elapsed time:', font=Fonts['button']).grid(row=1, column=0, sticky='e')
        labPid = tk.StringVar()
        tk.Label(frame, textvariable=labPid).grid(row=0, column=1, sticky='w')
        labElap = tk.StringVar()
        tk.Label(frame, textvariable=labElap).grid(row=1, column=1, sticky='w')
        tree.grid(row=2, column=0, sticky='nsew', columnspan=2)
        detailsSection.grid(row=1, column=1, sticky='nsew', padx=5, pady=5)

        # Handle close event
        root.protocol("WM_DELETE_WINDOW", onClose)
	
        # First time load in GUI thread        
        self.process = pm.getProcessFromPid()
        self.process_childs = pm.getProcessGroup()
        # Launch Working thread for next updates
        self.process_check = True
        self.process_thread = threading.Thread(target=processWorker)
        self.process_thread.start()
            
        refreshInfo()
        centerWindows(root, refWindows=self.root)
        root.deiconify()
        root.mainloop()
        
    def switchRunsView(self):
        '''Change the way to visualize protocol runs.
        It can be in a list view or tree view'''
        self.graphView = not self.graphView
        show = self.treeHist
        hide = self.graphHist
        iconPath = ActionIcons[ACTION_TREE]
        if self.graphView:
            show = self.graphHist
            hide = self.treeHist
            iconPath = "columns.gif"
        self.updateRunHistory(self.lastDisplayGroup)
        self.runButtonsDict[ACTION_TREE].setImage(iconPath)
        show.grid(row=0, column=0, sticky='nsew')
        hide.grid_remove()
    
    def launchDepsTree(self):
        runsDict = self.project.getRunsDependencies()
        from protlib_gui_figure import showDependencyTree
        #from protlib_gui_graph import showDependencyTree
        showDependencyTree(runsDict)
        
        
    def launchRunStepsTree(self, run):
        ''' Launch another windows with the steps 
        inserted for that run '''
        import pickle
        runName = getExtendedRunName(run)
        root = tk.Toplevel()
        root.withdraw()
        #root.title()
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)

        frame = tk.Frame(root)
        frame.grid(row=0, column=0, sticky='nsew')
        frame.columnconfigure(0, weight=1, minsize=400)
        frame.rowconfigure(1, weight=1)
        
        columns = ('Time')
        tree = XmippTree(frame, columns=columns)

        for c, name, width in [('#0', 'Step', 200), ('Time', 'Time', 100)]:
            tree.column(c, width=width)
            tree.heading(c, text=name)
            
        tree.grid(row=1, column=0, sticky='nsew')
        
        def insertSteps(e=None):
            tree.clear()
            r = tree.insert('', 'end', text=runName, image=self.getStateImage(run['run_state']))
            tree.item(r, open=True)
            steps = self.project.projectDb.getRunSteps(run)
            total = None
            
            for s in steps:    
                t = "%s (id = %d)" % (s['command'], s['step_id'])
                stepIcon = 'step_finished'
                
                if s['finish'] is None:
                    stepIcon = 'step'
                    timeValue = ''
                else:
                    init = datetime_fromdb(s['init']) 
                    end = datetime_fromdb(s['finish'])
                    delta = end - init
                    if total is None:
                        total = delta
                    else:
                        total += delta 
                    timeValue = pretty_delta(delta)
                    
                item = tree.insert(r, 'end', text=t, image=self.getImage(stepIcon), values=(timeValue,))
                p = tree.insert(item, 'end', text=' parameters', image=self.getImage('md_view'))
                vf = tree.insert(item, 'end', text=' verify_files', image=self.getImage('copy'))
                tree.insert(item, 'end', text='parent = %d' % s['parent_step_id'])
                tree.insert(item, 'end', text='   init = %s' % str(s['init']))
                tree.insert(item, 'end', text=' finish = %s' % str(s['finish']))
                
                params = pickle.loads(str(s['parameters']))
                for k, v in params.iteritems():
                    tree.insert(p, 'end', text=" %s = %s" % (k, v))
                
                files = pickle.loads(str(s['verifyFiles']))        
                for f in files:
                    tree.insert(vf, 'end', text=f)
            
            tree.item(r, values=(pretty_delta(total),))
            
        btn = XmippButton(frame, "Refresh", 'refresh.gif', command=insertSteps, tooltip='Refresh   F5')
        btn.grid(row=0, column=0, padx=(0, 5), sticky='nw')
        root.bind("<F5>", insertSteps)
        insertSteps()
        centerWindows(root, refWindows=self.root)
        root.deiconify()
        root.mainloop()       
        
    def launchProtocolGUI(self, run, visualizeMode=False):
        run['group_name'] = self.lastDisplayGroup
        top = tk.Toplevel()
        gui = ProtocolGUI()
        gui.createGUI(self.project, run, top, 
                      lambda: self.protocolSaveCallback(run), visualizeMode)
        gui.launchGUI()
        
    def protocolSaveCallback(self, run):
        #self.selectToolbarButton(run['group_name'], False)
        if self.lastDisplayGroup in [GROUP_ALL, run['group_name']]:
            self.historyRefreshRate = 1
            self.updateRunHistory(self.lastDisplayGroup)

    def updateRunHistory(self, protGroup, selectFirst=True, checkDead=False):
        #Cancel if there are pending refresh
        tree = self.treeHist
        graph = self.graphHist
        
        if self.historyRefresh:
            tree.after_cancel(self.historyRefresh)
            self.historyRefresh = None
        runName = None
        
        if not selectFirst and not self.graphView:
            runName = tree.item(tree.selection_first(), 'text')
        
        tree.clear()
        graph.clear()

        self.runs, stateList = self.project.getStateRunList(protGroup, checkDead)
        
        if len(stateList) > 0:
            if not self.graphView:
                # Update runs tree
                for name, state, stateStr, modified in stateList:
                    tree.insert('', 'end', text = name, tags=('normal'),
                                image=self.getStateImage(state),
                                values=(stateStr, modified)) 
                
                for c in tree.get_children(''):
                    if selectFirst or tree.item(c, 'text') == runName:
                        tree.selection_set(c)
                        self.updateRunSelection(tree.index(c), update_summary=checkDead)
                        break
            else:
                self.updateHistoryGraph()
            #Generate an automatic refresh after x ms, with x been 10, 100, 1000... increasing until 3600 (1h)
            self.historyRefreshRate = min(10*self.historyRefreshRate, 3600)
            self.historyRefresh = tree.after(self.historyRefreshRate*1000, self.updateRunHistory, protGroup, False)
        else:
            self.updateRunSelection(-1, update_summary=checkDead)

    #---------------- Functions related with Popup menu ----------------------   
    def unpostMenu(self, event=None):
        if self.lastMenu:
            self.lastMenu.unpost()
            self.lastMenu = None
        self.menuRun.unpost()
            
    def onRightClick(self, e=None):
        self.menuRun.delete(0, tk.END)
        aList = [(ACTION_EDIT, 'Edit     '),
                 (ACTION_COPY, 'Duplicate   '),
                 (ACTION_DELETE, 'Delete    '),
                 (None, None),
                 (ACTION_STOP, 'Stop'),
                 (ACTION_STEPS, 'Show steps')]
        
        def addMenuOption(action, label):
            if action is None: 
                self.menuRun.add_separator()
            else:
                imgName = ActionIcons[action].replace('.gif', '')
                self.menuRun.add_command(label=" "  + label, command=lambda: self.runButtonClick(action),
                                         image=self.getImage(imgName), compound=tk.LEFT)
        
        for action, label in aList:
            addMenuOption(action, label)
        self.menuRun.post(e.x_root, e.y_root)
                
    def clickToolbarButton(self, index, showMenu=True):
        key, btn, menu = self.ToolbarButtonsDict[index]
        self.unpostMenu()
        if menu:
            # Post menu
            x, y, w = btn.winfo_x(), btn.winfo_y(), btn.winfo_width()
            xroot, yroot = self.root.winfo_x() + btn.master.winfo_x(), self.root.winfo_y()+ btn.master.winfo_y()
            menu.post(xroot + x + w + 10, yroot + y)
            #if self.lastMenu:
            #    self.lastMenu.unpost()
            self.lastMenu = menu
            
    def cbDisplayGroupChanged(self, e=None):
        index = self.cbDisplayGroup.current()
        key, btn, menu = self.ToolbarButtonsDict[index]                
        self.selectDisplayGroup(key)
        
    def selectDisplayGroup(self, group=GROUP_ALL):
        '''Change the display group for runs list'''
        if not self.lastDisplayGroup or group != self.lastDisplayGroup:
            self.project.config.set('project', 'lastselected', group)
            self.project.writeConfig()
            self.updateRunHistory(group)
            self.lastDisplayGroup = group
            
    def updateRunSelection(self, index, update_summary=False):
        state = tk.NORMAL
        details = self.Frames['details']
        if index == -1:
            state = tk.DISABLED
            #Hide details
            details.grid_remove()
        else:
            state = tk.NORMAL
            run = self.lastRunSelected
            showButtons = False
            try:
                if not exists(run['script']):
                    return
                
                run_state = run['run_state']
                run_summary_cache = run['script'].replace('.py', '_summary.txt')
                update_summary = (update_summary or (run_state == SqliteDb.RUN_STARTED) or
                                  not exists(run_summary_cache))
                
                prot = self.project.getProtocolFromModule(run['script'])
                if update_summary:
                    comment = ""

                    if not prot.Comment == None:
                        if prot.Comment.find("Describe your run here...") < 0:
                            comment = prot.Comment
                            comment = comment.strip(' \t\n\r')
                            if len(comment) > 128:
                                comment = comment[0:128]+"..."
                    if exists(prot.WorkingDir):
                        summary = '\n'.join(prot.summary())
                        showButtons = True
                        wd = "[%s]" % prot.WorkingDir # If exists, create a link to open folder
                    else:
                        wd = prot.WorkingDir
                        summary = "This protocol run has not been executed yet"
                    
                    labels = '<Run>: ' + getExtendedRunName(run) + \
                              '\n<Created>: ' + run['init'] + '   <Modified>: ' + run['last_modified'] + \
                              '\n<Script>: ' + run['script'] + '\n<Directory>: ' + wd + \
                              '\n<Comment>: ' + comment + '\n\n<Summary>:\n' + summary
                    f = open(run_summary_cache, 'w')
                    f.write(labels)
                else:
                    
                    if exists(prot.WorkingDir):
                        showButtons = True
                    f = open(run_summary_cache)
                    labels = f.read()
                f.close()
                
            except Exception, e:
                labels = 'Error creating protocol: <%s>' % str(e)
                import traceback
                labels += '\n' + traceback.format_exc()
            self.detailsText.clear()
            self.detailsText.addText(labels)
            details.displayButtons(showButtons)
            details.grid()
        for btn in self.runButtonsDict.values():
            btn.config(state=state)
            
    def selectTreeRun(self, e=None):
        item = self.treeHist.selection_first()
        index = -1
        if item:
            index = self.treeHist.index(item)
            self.lastRunSelected = self.runs[index]
        self.updateRunSelection(index)
        
    def selectGraphRun(self, nodeName):
        if not nodeName.startswith('Project'):
            protName, runName = nodeName.split('\n')
            index = -1
            for i, run in enumerate(self.runs):
                if run['protocol_name'] == protName and run['run_name'] == runName:
                    index = i
                    self.lastRunSelected = run
                    break
            self.updateRunSelection(index)   
                        
    def doubleClickGraph(self, nodeName):
        self.runButtonClick(ACTION_DEFAULT)
        
    def rightClickGraph(self, event):
        nodeName = event.text
        if not nodeName.startswith('Project'):
            self.selectGraphRun(nodeName)
            self.onRightClick(event)
            
    def getLastRunDict(self):
        if self.lastRunSelected:
            from protlib_sql import getRunDict
            return getRunDict(self.lastRunSelected)
        return None
    
    def runButtonClick(self, event=None):
        #FlashMessage(self.root, 'Opening...', delay=1)
        run = self.getLastRunDict()
        if run:
            state = run['run_state']
            if event == ACTION_DEFAULT:
                if state == SqliteDb.RUN_STARTED:
                    self.launchRunJobMonitorGUI(run)
                else:
                    self.launchProtocolGUI(run)
            elif event == ACTION_EDIT:
                self.launchProtocolGUI(run)
            elif event == ACTION_COPY:
                self.launchProtocolGUI(self.project.copyProtocol(run['protocol_name'], run['script']))
            elif event == ACTION_DELETE:
                if state in [SqliteDb.RUN_STARTED, SqliteDb.RUN_LAUNCHED]:
                    showWarning("Delete warning", "This RUN is LAUNCHED or RUNNING, you need to stop it before DELETE", parent=self.root) 
                else:
                    runsDict = self.project.getRunsDependencies()
                    deps = ['<%s>' % d for d in runsDict[getExtendedRunName(run)].deps]
                    if len(deps):
                        showWarning("Delete warning", "This RUN is referenced from: \n" + '  \n'.join(deps), parent=self.root) 
                    else:
                        if askYesNo("Confirm DELETE", "<ALL DATA> related to this <protocol run> will be <DELETED>. \nDo you really want to continue?", self.root):
                            error = self.project.deleteRun(run)
                            if error is None:
                                self.updateRunHistory(self.lastDisplayGroup)
                            else: 
                                showError("Error deleting RUN", error, parent=self.root)
            elif event == ACTION_STEPS:
                self.launchRunStepsTree(run)
            elif event == ACTION_TREE:
                self.switchRunsView()
                #self.launchDepsTree()
            elif event == ACTION_REFRESH:
                self.historyRefreshRate = 1
                self.updateRunHistory(self.lastDisplayGroup, False, True)
            elif event == ACTION_STOP:
                self.stopRun()
        
    def createToolbarFrame(self, parent):
        #Configure toolbar frame
        toolbar = tk.Frame(parent, bd=2, relief=tk.RIDGE)
        
        self.Frames['toolbar'] = toolbar
        self.logo = getXmippImage('xmipp_logo.gif')
        label = tk.Label(toolbar, image=self.logo)
        label.pack()
        
        #Create toolbar buttons
        #i = 1
        section = None
        self.ToolbarButtonsDict[0] = (GROUP_ALL, None, None)
        index = 1
        for k, v in sections:
            section = ProjectButtonMenu(toolbar, k)
            section.pack(fill=tk.X, pady=5)
            for o in v:
                text = o[0]
                key = "%s_%s" % (k, text)
                opts = o[1:]
                def btn_command(index): 
                    def new_command(): 
                        self.clickToolbarButton(index)
                    return new_command 
                btn = section.addButton(text, command=btn_command(index))
                menu = self.createToolbarMenu(section, opts)
                self.ToolbarButtonsDict[index] = (key, btn, menu)
                index = index + 1
        #section.addButton(GROUP_XMIPP, command=self.launchProgramsGUI)
        #self.ToolbarButtonsDict[index] = (GROUP_XMIPP, None, None)

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
        
        history.addWidget(tk.Label, text="Filter")
        import ttk
        values = [GROUP_ALL]
        for k, v in sections:
            for o in v:
                values.append(o[0])
        self.cbDisplayGroup = history.addWidget(ttk.Combobox, values=values, state='readonly')
        self.cbDisplayGroup.bind('<<ComboboxSelected>>', self.cbDisplayGroupChanged)
        #cb.set(GROUP_ALL)
        
        aList = [(ACTION_EDIT, 'Edit    Alt-E/(Enter)'), 
                (ACTION_COPY, 'Duplicate   Alt-D/Ctrl-Enter'),
                (ACTION_DELETE, 'Delete    Del'),
                (ACTION_TREE, 'Show runs dependency tree'),
                (ACTION_REFRESH, 'Refresh   F5')]
        def setupButton(k, v, t):
            btn =  history.addButton(k, v, command=lambda:self.runButtonClick(k), bg=HighlightBgColor)
            ToolTip(btn, t, 500)
            self.runButtonsDict[k] = btn
        for k, v in aList:
            setupButton(k, ActionIcons[k], v)
        
        self.createHistoryTree(history.frameContent)  # set self.treeHist
        self.createHistoryGraph(history.frameContent) # set self.graphHist
        
        show = self.treeHist
        if self.graphView:
            show = self.graphHist
        show.grid(row=0, column=0, sticky='nsew')
        
        history.frameContent.columnconfigure(0, weight=1)
        history.frameContent.rowconfigure(0, weight=1)
        return history
    
    def createHistoryTree(self, parent):
        columns = ('State', 'Modified')
        tree = XmippTree(parent, columns=columns)
        for c in columns:
            tree.column(c, anchor='e')
            tree.heading(c, text=c) 
        tree.column('#0', width=300)
        tree.heading('#0', text='Run')
        tree.bind('<<TreeviewSelect>>', self.selectTreeRun)
        tree.bind('<Double-1>', lambda e:self.runButtonClick(ACTION_DEFAULT))
        tree.bind("<Button-3>", self.onRightClick) #right button on win and linux
        tree.bind("<Button-2>", self.onRightClick) #right button on mac
        tree.tag_configure('normal', font=Fonts['normal'])
        self.treeHist = tree
        
    def createHistoryGraph(self, parent):
        from protlib_gui_graph import Canvas
        canvas = Canvas(parent)
        canvas.onClickCallback = self.selectGraphRun
        canvas.onDoubleClickCallback = self.doubleClickGraph
        canvas.onRightClickCallback = self.rightClickGraph
        self.graphHist = canvas
        self.updateHistoryGraph()
        
    def updateHistoryGraph(self):
        runsDict = self.project.getRunsDependencies()
        #from protlib_gui_figure import showDependencyTree
        from protlib_gui_graph import showDependencyTree
        showDependencyTree(self.graphHist, runsDict, self.project.getName())
        
    def createDetailsFrame(self, parent):
        details = ProjectSection(parent, 'Details')
        self.Frames['details'] = details
        #Create RUN details
        details.addButton("Analyze results", command=self.visualizeRun, underline=0)
        details.addButton("Output files", command=self.showOutput, underline=0)
        content = details.frameContent
        content.config(bg=BgColor, bd=1, relief=tk.RIDGE)
        content.grid_configure(pady=(5, 0))
        
        self.detailsText = TaggedText(content, height=15, width=70)
        self.detailsText.frame.pack(fill=tk.BOTH)
        return details

    def createGUI(self, root=None):
        if not root:
            root = tk.Tk()
        self.root = root
        root.withdraw() # Hide the windows for centering
        projectName = self.project.getName()
        hostName = getHostname()        
        self.root.title("Xmipp Protocols   %(projectName)s  on  %(hostName)s" % locals())
        self.initVariables()        
        self.createMainMenu()

        #Create main frame that will contain all other sections
        #Configure min size and expanding behaviour
        root.minsize(750, 550)
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
        lastGroup = GROUP_ALL
        if self.project.config.has_option('project', 'lastselected'):
            lastGroup = self.project.config.get('project', 'lastselected')
        self.selectDisplayGroup(lastGroup)
        for k, v in self.ToolbarButtonsDict.iteritems():
            if lastGroup == v[0]:
                self.cbDisplayGroup.current(k)
                break
    
        self.addBindings()
                
    def launchGUI(self, center=True):
        if center:
            centerWindows(self.root)
        self.root.deiconify()
        self.root.mainloop()
       
    def close(self, event=""):
        self.root.destroy()

    def showOutput(self, event=''):
        prot = self.project.getProtocolFromModule(self.lastRunSelected['script'])
        title = "Output Console - %s" % self.lastRunSelected['script']
        filelist = ["%s%s" % (prot.LogPrefix, ext) for ext in ['.log', '.out', '.err']]
        showTextfileViewer(title, filelist, self.root)
        
    def visualizeRun(self, event=''):
        run = self.getLastRunDict()
        self.launchProtocolGUI(run, True)

from protlib_xmipp import XmippScript, greenStr, cyanStr, redStr

class ScriptProtocols(XmippScript):
    def __init__(self):
        XmippScript.__init__(self, True)
        
    def defineParams(self):
        self.addUsageLine("Create Xmipp project on this folder.");
        self.addUsageLine("Or manage an existing project");
        ## params
        self.addParamsLine("[ --list  ]                         : List project runs");
        self.addParamsLine("   alias -l;"); 
        self.addParamsLine("[ --details <runName>  ]            : Print details about some run");
        self.addParamsLine("   alias -e;"); 
        self.addParamsLine("[ --new_run <protName>  ]           : Create a new run of a protocol");
        self.addParamsLine("   alias -n;");
        self.addParamsLine("[ --copy_run <runName>  ]           : Create a copy of an existing run");
        self.addParamsLine("   alias -p;");
        self.addParamsLine("[ --start_run <runName> <mode=Resume> ] : Start(or launch to queue) an existing run");
        self.addParamsLine("   alias -s;");        
        self.addParamsLine("[ --delete_run <protName>  ]        : Delete an existing run");
        self.addParamsLine("   alias -d;");
        self.addParamsLine("[ --clean  ]                        : Clean ALL project data");
        self.addParamsLine("   alias -c;"); 
        self.addParamsLine("[ --update_project]                 : Update project db with new protocosl from 3.0 to 3.1")
        
    def confirm(self, msg, defaultYes=True):
        if defaultYes:
            suffix = ' [Y/n] '
        else:
            suffix = ' [y/N] '
        answer = raw_input(msg + cyanStr(suffix))
        if not answer:
            return defaultYes
        return answer.lower() == 'y'
    
    def printProjectName(self):
        '''Print project name '''
        print cyanStr("PROJECT: "), greenStr(os.path.basename(self.proj_dir))
        
    def loadProjectFromCli(self):
        '''Load project to work from command line'''
        self.printProjectName()
        self.launch = False
        self.project.load()

    def run(self):
        self.proj_dir = os.getcwd()
        self.project = XmippProject(self.proj_dir)
        self.launch = True

        if self.checkParam('--list'):
            self.loadProjectFromCli()
            runs, stateList = self.project.getStateRunList(checkDead=True)
            if len(stateList) > 0:
                print cyanStr("List of runs:")
                stateList.insert(0, ('RUN', 0, 'STATE', 'MODIFIED'))
                for name, state, stateStr, modified in stateList:
                    print "%s | %s | %s " % (name.rjust(35), stateStr.rjust(15), modified)
            else:
                print cyanStr("No runs found")
                
        elif self.checkParam('--details'):
            self.loadProjectFromCli()
            runName = self.getParam('--details')
            try:
                protocol = self.project.getProtocolFromRunName(runName)
                print cyanStr("Details of run "), greenStr(runName)
                summaryLines = protocol.summary()
                for line in summaryLines:
                    for c in '<>[]':
                        line = line.replace(c, '')
                    print "   ", line
            except Exception, e:
                print redStr("Run %s not found" % runName)
                
        elif self.checkParam('--new_run'):
            protName = self.getParam('--new_run')
            self.loadProjectFromCli()
            run = self.project.newProtocol(protName) #Create new run for this protocol
            self.project.projectDb.insertRun(run) # Register the run in Db
            ProtocolParser(run['source']).save(run['script']) # Save the .py file
            print "Run %s was created" % greenStr(getExtendedRunName(run))
            print "Edit file %s and launch it" % greenStr(run['script'])
            
        elif self.checkParam('--copy_run'):
            extRunName = self.getParam('--copy_run')
            self.loadProjectFromCli()
            protName, runName = splitExtendedRunName(extRunName)
            script = getScriptFromRunName(extRunName)
            run = self.project.copyProtocol(protName, script) #Copy parameter from existing run
            self.project.projectDb.insertRun(run) # Register the run in Db
            ProtocolParser(script).save(run['script']) # Save the .py file
            print "Run %s was created from %s" % (greenStr(getExtendedRunName(run)), greenStr(extRunName))
            print "Edit file %s and launch it" % greenStr(run['script'])
            
        elif self.checkParam('--delete_run'):
            extRunName = self.getParam('--delete_run')
            self.loadProjectFromCli()
            protName, runName = splitExtendedRunName(extRunName)
            self.project.deleteRunByName(protName, runName)
            print "Run deleted"      
            
        elif self.checkParam('--start_run'):
            self.loadProjectFromCli()
            extRunName = self.getParam('--start_run')
            print 'extRunName', extRunName
            mode = self.getParam('--start_run', 1) # read the mode argument
            protName, runName = splitExtendedRunName(extRunName)
            print 'protName, runName', protName, runName
            pe = ProtocolExecutor(protName, self.project, runName=runName)
            pe.setValue('Behavior', mode)
            pe.runProtocol()
            
        elif self.checkParam('--clean'):
            self.printProjectName()
            self.launch = self.confirm('%s, are you sure to clean?' % cyanStr('ALL RESULTS will be DELETED'), False)
            if self.launch:
                self.project.clean()
            else:
                print cyanStr("CLEAN aborted.")
                 
        elif self.checkParam('--update_project'):
            self.loadProjectFromCli()
            self.project.updateProtocolTables()
        else: #lauch project     
            if not self.project.exists():    
                print 'You are in directory: %s' % greenStr(self.proj_dir)
                self.launch = self.confirm('Do you want to %s in this folder?' % cyanStr('CREATE a NEW PROJECT'))
                if (self.launch):
                    self.project.create()
                    print "Created new project on directory: %s" % greenStr(self.proj_dir)
                else:
                    print cyanStr("PROJECT CREATION aborted")
            else:
                self.project.load()
        
        if self.launch:
            try:
                self.root = tk.Tk()
                self.root.withdraw()
                gui = XmippProjectGUI(self.project)
                gui.createGUI(self.root)
                gui.launchGUI()
            except Exception, e:
                print cyanStr("Could not create Tk GUI, disabled")
                import traceback
                import sys
                traceback.print_exc(file=sys.stderr)

if __name__ == '__main__':
    ScriptProtocols().tryRun()
