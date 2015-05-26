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
from pyworkflow.utils.properties import Message
import pyworkflow.gui as pwgui 
from pyworkflow.gui.widgets import HotButton
from pyworkflow.gui.text import TaggedText
from pyworkflow.gui.dialog import askString, askCheckString, askYesNo, showError
from pyworkflow.gui.browser import FileBrowserWindow

            
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
                     command=self.createNewProject)
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
            project = self.manager.loadProject(p.getName(), chdir=False)
            # Add creation time to project info
            p.cTime = project.getCreationTime()
            frame = self.createProjectLabel(parent, p, color=colors[i%2])
            frame.grid(row=r, column=0, padx=10, pady=5, sticky='new')
            r += 1
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
    
    def createNewProject(self, e=None):
        value =  askCheckString(Message.LABEL_CREATE_PROJECT, Message.TITLE_CREATE_PROJECT_NAME, self.root, "Change project location?", 30)
        if value is not None:
            (projName, askLocation) = value
            if not projName is None:
                def onCreate(projName, location = None):
                    self.manager.createProject(projName, location = location)
                    self.createProjectList(self.text)
                    self.openProject(projName)

                def onSelect(obj):
                    onCreate(projName, location=obj.getPath())

                if askLocation:
                    FileBrowserWindow("Browsing", self.windows, path=getHomePath(), onSelect=onSelect).show()
                else:
                    onCreate(projName)

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
