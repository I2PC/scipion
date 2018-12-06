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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This modules handles the System management
"""

import os

import pyworkflow as pw
import pyworkflow.utils as pwutils
from project import Project



class ProjectInfo(object):
    """Class to store some information about the project"""
    def __init__(self, projName, mTime):
        """At least it receives the Project Name and its modification time"""
        self.projName = projName
        self.mTime = mTime
        
    def getName(self):
        return self.projName
    
    def getModificationTime(self):
        return self.mTime
        
        
class Manager(object):
    """ Manage all projects of a given workspace.
     (i.e, listing projects, creating, deleting or querying info)
    """
    def __init__(self, workspace=None):
        """
        Params:
            workspace: path to where the user's workspace is. A subfolder
            named 'projects' is expected. If workspace is None, the workspace
            will be taken from the global configuration.
        """
        ws = workspace or pw.Config.SCIPION_USER_DATA
        self.PROJECTS = os.path.join(ws, 'projects')

    def getProjectPath(self, projectName):
        """Return the project path given the name"""
        return os.path.join(self.PROJECTS, projectName)
        
    def listProjects(self, sortByDate=True):
        """Return a list with all existing projects
        And some other project info
        If sortByData is True, recently modified projects will be first"""
        projList = []
        pw.utils.makePath(self.PROJECTS)
        for f in os.listdir(self.PROJECTS):
            p = self.getProjectPath(f)
            if os.path.isdir(p):
                stat = os.stat(p)
                projList.append(ProjectInfo(f, stat.st_mtime))
                
        if sortByDate:
            projList.sort(key=lambda k: k.mTime, reverse=True)
        return projList
    
    def createProject(self, projectName, runsView=1, 
                      hostsConf=None, protocolsConf=None, location=None):
        """Create a new project.
        confs dict can contains customs .conf files 
        for: menus, protocols, or hosts
        """
        # If location is not None create project on it (if exists)
        if location is None:
            projectPath = self.getProjectPath(projectName)
        else:
            projectPath = os.path.join(location, projectName)

        # JMRT: Right now the project.create function change the current
        # working dir (cwd) to the project location, let's store the cwd
        # and restored after the creation
        cwd = os.getcwd()
        project = Project(projectPath)
        project.create(runsView=runsView, 
                       hostsConf=hostsConf, 
                       protocolsConf=protocolsConf)
        # If location is not the default one create a symlink on self.PROJECTS directory
        if projectPath != self.getProjectPath(projectName):
            # JMRT: Let's create the link to the absolute path, since relative
            # can be broken in systems with different mount points
            pw.utils.createAbsLink(os.path.abspath(projectPath),
                                   self.getProjectPath(projectName))

        os.chdir(cwd)  # Retore cwd before project creation

        return project

    def importProject(self, fromLocation, copyFiles=True, projectName=None, searchLocation=None):
        """Import a project that is somewhere else in the FS
        Folder can be copied (default) or linked
        Optionally a name can be specified, otherwise name will match location folder name
        Project will always be created in the default project folder.
        """
        # If projectName is None...
        if projectName is None:
            # use same name as the import location
            projectName = os.path.basename(fromLocation)

        projectPath = self.getProjectPath(projectName)

        # If need to copyFiles
        if copyFiles:
            # Copy the whole folder
            pw.utils.copyTree(os.path.abspath(fromLocation), os.path.abspath(projectPath))

        else:
            # Link the folder
            pw.utils.createAbsLink(os.path.abspath(fromLocation), projectPath)

        project = self.loadProject(projectName)

        if searchLocation: project.fixLinks(searchLocation)

        return project

    def loadProject(self, projId, **kwargs):
        """ Retrieve a project object, given its id. """
        project = Project(self.getProjectPath(projId))
        project.load(**kwargs)
        return project

    def deleteProject(self, projectName):
        pw.utils.cleanPath(self.getProjectPath(projectName))

    def renameProject(self, oldName, newName):
        os.rename(self.getProjectPath(oldName), self.getProjectPath(newName))

    def hasProject(self, projectName):
        """Return True if exists a project with projectName"""
        for projInfo in self.listProjects():
            if projectName == projInfo.projName:
                return True
        return False
