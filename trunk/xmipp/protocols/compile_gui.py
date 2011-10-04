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
from Tkinter import *
import tkFont as font
import ttk  

from protlib_gui_ext import FilePollTextArea, centerWindows, \
                            MyButton, registerCommonFonts, showInfo, showError
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
class OptionsTab(Frame):
    def __init__(self, master, text, **opts):
        Frame.__init__(self, master)
        self.config(**opts)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=3)
        self.lastRow = 0
        self.options = []
        self.optionsDict = {}
        style = ttk.Style()
        style.configure("Name.TLabel", foreground="red")
        self.normal = font.Font(family='Helvetica', size=10)
        self.bold = font.Font(family='Helvetica', size=10, weight='bold')
        
    def addOption(self, name, comment, default='', cond=None, wiz=None, browse=False):
        r = self.lastRow
        var = StringVar()
        var.set(default)
        ttk.Label(self, text=comment, font=self.normal).\
        grid(column=0, row=r, padx=5, pady=5, sticky=E)
        ttk.Label(self, text=name, font=self.bold).\
        grid(column=1, row=r, padx=5, pady=5, sticky=E)
        if default in ['yes', 'no']: # Boolean option
            w = ttk.Checkbutton(self, textvariable=var, variable=var,
                            onvalue='yes', offvalue='no',
                            command=lambda:self.checked(name, var))
            w.grid(column=2, row=r, padx=5, pady=5, sticky=W)
        else:
            w = ttk.Entry(self, width=20, textvariable=var)
            
            if cond:
                if self.optionsDict[cond].get() == "no":
                    w['state'] = 'disabled'
            w.grid(column=2, row=r, sticky=(W, E), padx=5, pady=5)
            
            if browse:
                btn = MyButton(self, 'Browse', 'folderopen.gif',
                               command=lambda: browseDir(var, self))
                btn.grid(column=3, row=r)
                
            if wiz:
                btn = MyButton(self, 'Find', 'wizard.gif',
                               command=lambda: detectDir(self, wiz))
                btn.grid(column=4, row=r)
                
        self.options.append((name, default, var, w, cond))
        self.optionsDict[name] = var
        self.lastRow += 1
        
    def setValue(self, name, value):
        self.optionsDict[name].set(value)
        
    def getValue(self, name):
        return self.optionsDict[name].get()
        
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
        ttk.Separator(self, orient=HORIZONTAL).\
        grid(column=0, row=self.lastRow, columnspan=3,
             sticky=(W, E), padx=10, pady=5)
        self.lastRow += 1
        
    def getConfigOptions(self):
        optStr = ""
        for n, d, v, w, cond in self.options:
            optStr += ' %s="%s"' % (n, v.get())
        return optStr
        
class ConfigNotebook(ttk.Notebook):
    def __init__(self, master, OUTPUT, options, runFunc, stopCompileFunc):
        ttk.Notebook.__init__(self, master)
        self.tabs = {}
        self.master = master
        self.options = options
        self.run = runFunc
        self.stop = stopCompileFunc
        self.OUTPUT = OUTPUT
        
    def addTab(self, text):
        tab = OptionsTab(self, text)
        self.add(tab, text=text)
        self.tabs[text] = tab
        return tab
    
    def getConfigOptions(self):
        '''return string with configs settings '''
        return ' '.join([t.getConfigOptions() for t in self.tabs.values()])
    
    def setValue(self, tab, option, value):
        self.tabs[tab].setValue(option, value)
        
    def getValue(self, tab, option):
        return self.tabs[tab].getValue(option)
    
    def notifyRun(self, process):
        self.proc = process
        self.text.fillTextArea()
        self.text.doRefresh(1)
        self.checkProcess()    
        self.btn.config(command=self.notifyStopCompile, text="Stop")
        self.select(self.index('end') - 1)
        
    def notifyStopCompile(self, msg, isError=False):
        runFunc = self.run
        self.btn.config(command=lambda:launchRun(self, runFunc), text="Compile")
        
    def createConfigTab(self):
        tab = self.addTab("Configure")
        self.text = FilePollTextArea(tab, self.OUTPUT, 20, 80, colorOn=False)
        self.text.goEnd()
        self.text.grid(column=0, row=0, sticky=(N, S, E, W), padx=10, pady=10)

    def checkProcess(self):
        if self.proc.poll() is not None:
            self.text.stopRefresh()
            self.text.fillTextArea(goEnd=True)
            self.progressVar.set(0)
            if self.proc.returncode != 0:
                showError("Errors", "Errors on Xmipp compilation, see <%s> for more details" % self.OUTPUT, self)
            else:
                infoMsg = "<Xmipp> has been successfully compiled     \n"
                if self.options.hasOption('install'):
                    infoMsg += "<INSTALLATION FINISHED!!!>\n"
                    infoMsg += "Include file <.xmipp.bashrc> or <.xmipp.csh> to your startup shell file"
                showInfo("Compilation FINISHED", infoMsg, self)
                self.master.destroy()
        else:
            self.master.after(3000, self.checkProcess)
            if os.path.exists(self.OUTPUT):
                lines = 0
                for line in open(self.OUTPUT):
                    lines += 1
                self.progressVar.set(lines)
                
    def createPanels(self, runFunc):
        root = self.master
        #left panel
        leftFrame = ttk.Frame(root)
        imgPath = getXmippPath('resources', 'xmipp_logo.gif')
        self.img = PhotoImage(file=imgPath)
        label = ttk.Label(leftFrame, image=self.img)
        label.grid(column=0, row=0)
        leftFrame.grid(column=0, row=0, sticky=(N, S), padx=5, pady=5, rowspan=2)
        self.grid(column=1, row=0, sticky=(N, W, E, S), padx=5, pady=5)
        #bottom panel
        panel = ttk.Frame(root)
        panel.grid(column=1, row=1, padx=5, pady=5, sticky=[W, E])
        self.progressVar = IntVar()
        self.progressVar.set(0)
        progress = ttk.Progressbar(panel, orient=HORIZONTAL, length=300, mode='determinate', variable=self.progressVar, maximum="700")
        progress.pack(side=LEFT, padx=(0, 30))
        
        registerCommonFonts()
        self.btn = MyButton(panel, text='Compile')
        self.btn.pack(side=RIGHT, padx=(15, 0))
        procVar = StringVar()
        procVar.set(self.options.getNumberOfCpu())
        self.procVar = procVar
        procEntry = ttk.Entry(panel, width="5", textvariable=procVar)
        procEntry.pack(side=RIGHT, padx=5)
        Label(panel, text="Processors").pack(side=RIGHT, padx=5)
        btnFunc = lambda:launchRun(self, runFunc)
        self.btn.config(command=btnFunc)
        self.btn.bind('<Return>', func=lambda e:btnFunc())
        self.btn.focus_set()

def launchRun(nb, runFunc):
    nb.options.setNumberOfCpu(nb.procVar.get())
    runFunc(nb)
   
def createGUINotebook(OUTPUT, options, addTabsFunc, runFunc, stopCompileFunc):
    root = Tk()
    root.withdraw()
    root.title("Xmipp Install")
    root.minsize(width=600, height=350)
    nb = ConfigNotebook(root, OUTPUT, options, runFunc, stopCompileFunc)
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
    
