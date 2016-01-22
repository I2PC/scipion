# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'xmipp@cnb.csic.es'
# *
# **************************************************************************

import os
import Tkinter as tk
import tkFont

import pyworkflow as pw
from pyworkflow.utils.utils import prettyDate, prettyTime
from pyworkflow.utils.path import getHomePath
from pyworkflow.manager import Manager
import pyworkflow.gui as pwgui
from pyworkflow.gui.text import TaggedText
from pyworkflow.gui.dialog import askString, askYesNo, showError

from pyworkflow.gui import Message, Window, cfgEntryBgColor
from pyworkflow.gui.browser import FileBrowserWindow
from pyworkflow.gui.widgets import IconButton, HotButton, Button
from pyworkflow.utils.properties import Icon

            
class ProjectsView(tk.Frame):    
    def __init__(self, parent, windows, **args): 
        tk.Frame.__init__(self, parent, bg='white', **args)
        self.windows = windows
        self.manager = windows.manager
        self.root = windows.root
        
        #tkFont.Font(size=12, family='verdana', weight='bold')
        bigSize = pwgui.cfgFontSize + 2
        smallSize = pwgui.cfgFontSize - 2
        fontName = pwgui.cfgFontName
        
        self.projNameFont = tkFont.Font(size=bigSize, family=fontName, weight='bold')
        self.projDateFont = tkFont.Font(size=smallSize, family=fontName)
        self.projDelFont = tkFont.Font(size=smallSize, family=fontName, weight='bold')
        self.manager = Manager()
        btn = HotButton(self, text=Message.LABEL_CREATE_PROJECT, font=self.projNameFont, 
                     command=self._onCreateProject)
        btn.grid(row=0, column=0, sticky='nw', padx=10, pady=10)
        
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        text = TaggedText(self, width=40, height=15, bd=0, bg='white')
        text.grid(row=1, column=0, sticky='news')
      
        self.createProjectList(text)
        text.setReadOnly(True)
        self.text = text
    
    def createProjectList(self, text):
        """Load the list of projects"""
        r = 0
        text.setReadOnly(False)
        text.clear()
        parent = tk.Frame(text, bg='white')    
        parent.columnconfigure(0, weight=1)
        colors = ['white', '#EAEBFF']
        for i, p in enumerate(self.manager.listProjects()):
            try:
                project = self.manager.loadProject(p.getName(), chdir=False)
                # Add creation time to project info
                p.cTime = project.getCreationTime()
                frame = self.createProjectLabel(parent, p, color=colors[i%2])
                frame.grid(row=r, column=0, padx=10, pady=5, sticky='new')
                r += 1
            except Exception, ex:
                print "ERROR loading project: %s" % p.getName()
                print ex
        text.window_create(tk.INSERT, window=parent)
        text.bindWidget(parent)
        text.setReadOnly(True)
      
    def createProjectLabel(self, parent, projInfo, color):
        frame = tk.Frame(parent, bg=color)
        label = tk.Label(frame, text=projInfo.projName, anchor='nw', bg=color, 
                         justify=tk.LEFT, font=self.projNameFont, cursor='hand1', width=50)
        label.grid(row=0, column=0,  padx=2, pady=2, sticky='nw')
        label.bind('<Button-1>', lambda e: self.openProject(projInfo.projName))
        dateMsg = '%s%s    %s%s' % (Message.LABEL_MODIFIED, prettyDate(projInfo.mTime),
                                    Message.LABEL_CREATED, prettyTime(projInfo.cTime, time=False))
        dateLabel = tk.Label(frame, text=dateMsg, font=self.projDateFont, bg=color)
        dateLabel.grid(row=1, column=0, sticky='nw')
        delLabel = tk.Label(frame, text=Message.LABEL_DELETE_PROJECT, font=self.projDelFont, bg=color, cursor='hand1')
        delLabel.grid(row=1, column=1, padx=10)
        delLabel.bind('<Button-1>', lambda e: self.deleteProject(projInfo.projName))
        mvLabel = tk.Label(frame, text=Message.LABEL_RENAME_PROJECT, font=self.projDelFont, bg=color, cursor='hand1')
        mvLabel.grid(row=1, column=2)
        mvLabel.bind('<Button-1>', lambda e: self.renameProject(projInfo.projName))
        
        return frame
    
    def createNewProject(self, projName, projLocation):
        self.manager.createProject(projName, location=projLocation)
        self.createProjectList(self.text)
        self.openProject(projName)

    def _onCreateProject(self, e=None):
        projWindow = ProjectCreateWindow("Create project", self)
        projWindow.show()

    def openProject(self, projName):
        from subprocess import Popen
        script = pw.join('apps', 'pw_project.py')
        Popen([os.environ['SCIPION_PYTHON'], script, projName])
          
    def deleteProject(self, projName):
        if askYesNo(Message.TITLE_DELETE_PROJECT, 
                    "Project *%s*. " % projName + Message.MESSAGE_DELETE_PROJECT, self.root):
            self.manager.deleteProject(projName)
            self.createProjectList(self.text)

    def renameProject(self, projName):
        newName = askString("Rename project %s" % projName, "Enter new name:", self.root)
        if not newName or newName == projName:
            return
        if self.manager.hasProject(newName):
            showError("Rename cancelled",
                      "Project name already exists: %s" % newName, self.root)
            return
        self.manager.renameProject(projName, newName)
        self.createProjectList(self.text)


class ProjectCreateWindow(Window):
    """ Windows to create a project. """
    def __init__(self, title, parent=None, weight=True, minsize=(400, 150),
                 icon="scipion_bn.xbm", **args):
        """
         We assume the parent should be of ProjectsView
        """
        Window.__init__(self, title, parent.windows, weight=weight,
                        icon=icon, minsize=minsize, enableQueue=True)

        self.parent = parent
        self.projectsPath = self.parent.manager.PROJECTS
        self.projName = tk.StringVar()
        self.projName.set('')
        self.projLocation = tk.StringVar()
        self.projLocation.set(self.projectsPath)

        content = tk.Frame(self.root)
        content.columnconfigure(0, weight=1)
        content.columnconfigure(1, weight=1)
        content.config(bg='white')
        content.grid(row=0, column=0, sticky='news',
                       padx=5, pady=5)
        labelName = tk.Label(content, text=Message.LABEL_PROJECT, bg='white', bd=0)
        labelName.grid(row=0, column=0, sticky='nw', padx=5, pady=5)
        entryName = tk.Entry(content, bg=cfgEntryBgColor, width=20, textvariable=self.projName)
        entryName.grid(row=0, column=1, sticky='nw', padx=5, pady=5)

        labelCheck = tk.Label(content, text="Use default location", bg='white', bd=0)
        labelCheck.grid(row=1, column=0, sticky='nw', padx=5, pady=5)
        self.tkCheckVar = tk.IntVar()
        btnCheck = tk.Checkbutton(content, variable=self.tkCheckVar, bg='white', bd=0)
        btnCheck.grid(row=1, column=1, sticky='nw', padx=5, pady=5)

        self.browseFrame = tk.Frame(content, bg='white')
        #self.browseFrame.columnconfigure(1, weight=1)
        self.browseFrame.grid(row=2, column=0, padx=0, pady=0, columnspan=2, sticky='nw')
        self.entryBrowse = tk.Entry(self.browseFrame, bg=cfgEntryBgColor, width=40, textvariable=self.projLocation)
        self.entryBrowse.grid(row=0, column=0, sticky='nw', padx=5, pady=5)
        self.btnBrowse = IconButton(self.browseFrame, 'Browse', Icon.ACTION_BROWSE, command=self._browsePath)
        self.btnBrowse.grid(row=0, column=1, sticky='e', padx=5, pady=5)

        self.initial_focus = entryName
        self.tkCheckVar.trace('w', self._onVarChanged)
        btnCheck.select()

        btnFrame = tk.Frame(content)
        btnFrame.columnconfigure(0, weight=1)
        btnFrame.grid(row=3, column=0, sticky='sew', padx=5, pady=(0, 5), columnspan=2)
        btnFrame.config(bg='white')

        # Create buttons
        btnSelect = HotButton(btnFrame, 'Create', Icon.BUTTON_SELECT, command=self._select)
        btnSelect.grid(row=0, column=0, sticky='e', padx=5, pady=5)
        btnCancel = Button(btnFrame, 'Cancel', Icon.BUTTON_CANCEL, command=self.close)
        btnCancel.grid(row=0, column=1, sticky='e', padx=5, pady=5)


    def _onVarChanged(self, *args):
        if self.tkCheckVar.get() == 0:
            self.entryBrowse.config(state='normal')
            self.btnBrowse.config(state='normal')
        else:
            self.projLocation.set(self.projectsPath)
            self.entryBrowse.config(state='readonly')
            self.btnBrowse.config(state='disabled')

    def _browsePath(self, e=None):
        def onSelect(obj):
            self.projLocation.set(obj.getPath())

        v = self.projLocation.get().strip()
        path = None
        if v:
            v = os.path.dirname(v)
            if os.path.exists(v):
                path = v
        if not path:
            path = self.projectsPath

        browser = FileBrowserWindow("Browsing", self, path=path, onSelect=onSelect)
        browser.show()

    def _select(self):
        projName = self.projName.get().strip()
        projLocation = self.projLocation.get().strip()

        # Validate that project name is not empty
        if not projName:
            showError("Validation error", "Project name is empty", self.root)
        # Validate that project location is not empty
        elif not projLocation:
            showError("Validation error", "Project location is empty", self.root)
        # Validate that project location exists
        elif not os.path.exists(projLocation):
            showError("Validation error", "Project location does not exist", self.root)
        # Validate that project location is a directory
        elif not os.path.isdir(projLocation):
            showError("Validation error", "Project location is not a directory", self.root)
        # Validate that project path (location + name) ddoes not exists
        elif os.path.exists(os.path.join(projLocation, projName)):
            showError("Validation error", "Project path already exists", self.root)
        else:
            self.parent.createNewProject(projName, projLocation)
            self.close()
