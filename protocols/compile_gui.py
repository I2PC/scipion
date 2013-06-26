#!bin/xmipp_python
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
import tkFont as font
import ttk  

from protlib_gui_ext import centerWindows, XmippButton, registerCommonFonts, showInfo, showError, OutputText,\
    Fonts
from config_protocols import FontName, FontSize
from protlib_filesystem import getXmippPath
    
############### Helper functions ######################
def browseDir(var, parent):
    import tkFileDialog
    path = tkFileDialog.askdirectory(title="Choose directory", parent=parent)
    if len(path) > 0:
        var.set(os.path.abspath(path))
        
def detectDir(tab, detectFunc):
    errors = detectFunc(tab)
    if len(errors) > 0:
        showError("Errors", errors, tab)

################ Classes for a GUI based configuration ######################
class OptionsTab(tk.Frame):
    def __init__(self, master, text, **opts):
        tk.Frame.__init__(self, master)
        self.config(**opts)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=3)
        self.lastRow = 0
        self.options = []
        self.optionsDict = {}
	self.optionsGroup = {}
	self.optionsGroupPanels = {}
        style = ttk.Style()
        style.configure("Name.TLabel", foreground="red")
        self.normal = font.Font(family=FontName, size=FontSize)
        self.bold = font.Font(family=FontName, size=FontSize, weight='bold')
        
    def addOption(self, name, comment, default='', group=None, cond=None, wiz=None, browse=False):
        r = self.lastRow
        var = tk.StringVar()
        var.set(default)
	selfish = self
	if group != None:
	    selfish = self.getGroupPanel(group)
        l1 = ttk.Label(selfish, text=comment, font=self.normal).grid(column=0, row=r, padx=5, pady=5, sticky='e')
        l2 = ttk.Label(selfish, text=name, font=self.bold).grid(column=1, row=r, padx=5, pady=5, sticky='e')
        if default in ['yes', 'no']: # Boolean option
            w = ttk.Checkbutton(selfish, textvariable=var, variable=var,
                            onvalue='yes', offvalue='no',
                            command=lambda:self.checked(name, var))
            w.grid(column=2, row=r, padx=5, pady=5, sticky='w')
        else:
            w = ttk.Entry(selfish, width=20, textvariable=var)
            
            if cond:
                if self.optionsDict[cond].get() == "no":
                    w['state'] = 'disabled'
            w.grid(column=2, row=r, sticky='we', padx=5, pady=5)
            
            if browse:
                btn = XmippButton(selfish, 'Browse', 'folderopen.gif',
                               command=lambda: browseDir(var, self))
                btn.grid(column=3, row=r)
                
            if wiz:
                btn = XmippButton(selfish, 'Find', 'wizard.gif',
                               command=lambda: detectDir(self, wiz))
                btn.grid(column=4, row=r)
                
        self.options.append((name, default, var, w, cond))
        self.optionsDict[name] = var
        self.lastRow += 1

    def setValue(self, name, value):
        self.optionsDict[name].set(value)

    def setGroup(self, name, group):
        self.optionsGroup[name].set(group)

    def getValue(self, name, group=None):
        if group != None:
            return self.optionsDict[name].get()
        return self.optionsGroupPanels[group].optionsDict[name].get()

#    def getGroup(self, name):
#        return self.optionsGroup[name].get()

    def checked(self, name, var):
        value = var.get()
        if value == 'yes':
            state = 'normal'
        else:
            state = 'disabled'        
        for n, d, v, w, cond in self.options:
            if name == cond:
                w['state'] = state
                
    def addSeparator(self):
        ttk.Separator(self, orient=tk.HORIZONTAL).\
        grid(column=0, row=self.lastRow, columnspan=3,
             sticky='we', padx=10, pady=5)
        self.lastRow += 1

    def addGroupPanel(self, title):
        p = tk.LabelFrame(self, text=title)
	self.optionsGroupPanels[title] = p
	p.pack()

    def finishGroupPanel(self, title):
        p = self.getGroupPanel(title)
	p.pack(fill="both", expand="yes")

    def getGroupPanel(self, title):
        return self.optionsGroupPanels[title]

    def getConfigOptions(self):
        optStr = ""
        for n, d, v, w, cond in self.options:
            optStr += ' %s="%s"' % (n, v.get())
        return optStr
        
class ConfigNotebook(ttk.Notebook):
    def __init__(self, master, OUTPUT, options, runFunc):
        ttk.Notebook.__init__(self, master)
        self.tabs = {}
        self.master = master
        self.options = options
        self.run = runFunc
        self.OUTPUT = OUTPUT
        
    def addTab(self, text):
        tab = OptionsTab(self, text)
        self.add(tab, text=text)
        self.tabs[text] = tab
        return tab

    def getTab(self, text):
        return self.tabs[text]

    def getConfigOptions(self):
        '''return string with configs settings '''
        return ' '.join([t.getConfigOptions() for t in self.tabs.values()])
    
    def setValue(self, tab, option, value):
        self.tabs[tab].setValue(option, value)
        
    def getValue(self, tab, option, group=None):
        return self.tabs[tab].getValue(option, group)

    def notifyRun(self, process):
        self.proc = process
        self.text.readFile()
        self.text.doRefresh(1)
        self.checkProcess()    
        self.btn.config(command=lambda:self.stopCompile(), text="Stop")
        self.select(self.index('end') - 1)
        
    def stopCompile(self):
        self.master.after_cancel(self.checkRefresh)
        self.proc.terminate()
        runFunc = self.run
        self.btn.config(command=lambda:launchRun(self, runFunc), text="Compile")
        showInfo("Compilation STOPPED", "Compilation has been aborted.     ", self)
        
    def createConfigTab(self):
        tab = self.addTab("  Output  ")
        tab.columnconfigure(0, weight=1)
        tab.rowconfigure(0, weight=1)
        self.text = OutputText(tab, self.OUTPUT, width=80, height=20, colors=False)
        self.text.goEnd()
        self.text.frame.grid(column=0, row=0, sticky='nsew', padx=10, pady=10)

    def checkProcess(self):
        if self.proc.poll() is not None:
            self.text.stopRefresh()
            self.text.readFile()
            self.progressVar.set(0)
            if self.proc.returncode != 0:
                runFunc = self.run
                self.btn.config(command=lambda:launchRun(self, runFunc), text="Compile")
                showError("Errors", "Errors on <Xmipp> compilation, see <%s> for more details" % self.OUTPUT, self)
            else:
                infoMsg = "<Xmipp> has been successfully compiled     \n"
                if self.options.hasOption('install'):
                    infoMsg += "<INSTALLATION FINISHED!!!>\n"
                    infoMsg += "Include file <.xmipp.bashrc> or <.xmipp.csh> to your startup shell file"
                showInfo("Compilation FINISHED", infoMsg, self)
                self.master.destroy()
        else:
            self.checkRefresh = self.master.after(3000, self.checkProcess)
            if os.path.exists(self.OUTPUT):
                lines = 0
                for line in open(self.OUTPUT):
                    lines += 1
                self.progressVar.set(lines)
                
    def createPanels(self, runFunc):
        root = self.master
        root.columnconfigure(0, minsize=130, weight=1)
        registerCommonFonts()
        #left panel
        bgColor = 'white'
        leftFrame = tk.Frame(root, bg=bgColor)
        leftFrame.columnconfigure(0, weight=1)
        leftFrame.rowconfigure(0, minsize=100)
        imgPath = getXmippPath('resources', 'xmipp_logo.gif')
        self.img = tk.PhotoImage(file=imgPath)
        tk.Label(leftFrame, image=self.img, bg=bgColor).grid(column=0, row=0, sticky='we')
        tk.Label(leftFrame, text='Xmipp 3.0',  font=Fonts['button'], bg=bgColor).grid(column=0, row=1, sticky='we')
#       TODO: insert label extracting it from git repository
#        tk.Label(leftFrame, text='r12.4.3.11834', bg=bgColor).grid(column=0, row=2, sticky='we')
        leftFrame.grid(column=0, row=0, sticky='nsew', padx=5, pady=5, rowspan=2)
        self.grid(column=1, row=0, sticky='nsew', padx=5, pady=5)
        #bottom panel
	panel = ttk.Panedwindow(root)
        #panel = ttk.Frame(root)
        panel.grid(column=1, row=1, padx=5, pady=5, sticky='we')
        self.progressVar = tk.IntVar()
        self.progressVar.set(0)
        progress = ttk.Progressbar(panel, orient=tk.HORIZONTAL, length=300, mode='determinate', variable=self.progressVar, maximum="700")
        progress.pack(side=tk.LEFT, padx=(0, 30))
        
        XmippButton(panel, text='Close', command=lambda: self.master.destroy()).pack(side=tk.RIGHT, padx=(5, 0))
        self.btn = XmippButton(panel, text='Compile')
        self.btn.pack(side=tk.RIGHT, padx=(15, 0))
        procVar = tk.StringVar()
        procVar.set(self.options.getNumberOfCpu())
        self.procVar = procVar
        procEntry = ttk.Entry(panel, width="5", textvariable=procVar)
        procEntry.pack(side=tk.RIGHT, padx=5)
        tk.Label(panel, text="Processors").pack(side=tk.RIGHT, padx=5)
        btnFunc = lambda:launchRun(self, runFunc)
        self.btn.config(command=btnFunc)
        self.btn.bind('<Return>', func=lambda e:btnFunc())
        self.btn.focus_set()

def launchRun(nb, runFunc):
    nb.options.setNumberOfCpu(nb.procVar.get())
    runFunc(nb)
   
def createGUINotebook(OUTPUT, options, addTabsFunc, runFunc):
    from protlib_utils import getHostname
    root = tk.Tk()
    root.withdraw()
    root.title("Xmipp Install on " + getHostname())
    root.minsize(width=700, height=440)
    nb = ConfigNotebook(root, OUTPUT, options, runFunc)
    addTabsFunc(nb)
    nb.createConfigTab()
    root.columnconfigure(0, weight=1)
    root.columnconfigure(1, weight=5)
    root.rowconfigure(0, weight=1)
    nb.createPanels(runFunc)
    centerWindows(root)
    root.deiconify()
    root.mainloop()
    return nb
    
