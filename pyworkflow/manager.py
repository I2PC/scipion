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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This modules handles the System management
"""
import os
from os.path import abspath, join

from project import Project
import pyworkflow as pw
from pyworkflow.mapper import SqliteMapper
from pyworkflow.utils.path import cleanPath, makePath, getHomePath, existsPath
from pyworkflow.hosts import HostMapper
from pyworkflow.apps.config import writeHosts, getSettingsPath


SCIPION_PATH = 'Scipion'
PROJECTS_PATH = 'projects'


class ProjectInfo(object):
    """Class to store some information about the project"""
    def __init__(self, projName, mTime):
        """At least it receives the Project Name and its modification time"""
        self.projName = projName
        self.mTime = mTime
        
        
class Manager(object):
    """This class will handle the creation, modification
    and listing of projects."""
    def __init__(self):
        """For create a Project, the path is required"""
        self.path = join(getHomePath(), SCIPION_PATH, PROJECTS_PATH)
        makePath(self.path)
        settingsPath = getSettingsPath() 
        if not existsPath(settingsPath):
            writeHosts(settingsPath)
        
    def getProjectPath(self, projectName):
        """Return the project path given the name"""
        return join(self.path, projectName)
        
    def listProjects(self, sortByDate=True):
        """Return a list with all existing projects
        And some other project info
        If sortByData is True, recently modified projects will be first"""
        projList = []
        for f in os.listdir(self.path):
            p = self.getProjectPath(f)
            if os.path.isdir(p):
                stat = os.stat(p)
                projList.append(ProjectInfo(f, stat.st_mtime))
                
        if sortByDate:
            projList.sort(key=lambda k: k.mTime, reverse=True)
        return projList
    
    def createProject(self, projectName):
        """Create a new project """
        proj = Project(self.getProjectPath(projectName))
        hosts = HostMapper(getSettingsPath()).selectAll()
        proj.create(hosts)
        return proj
    
    def loadProject(self, projId):
        """ Retrieve a project object, given its id. """
        proj = Project(self.getProjectPath(projId))
        proj.load()
        return proj
        
    def deleteProject(self, projectName):
        cleanPath(self.getProjectPath(projectName))
        
    def hasProject(self, projectName):
        """Return True if exists a project with projectName"""
        for projInfo in self.listProjects():
            if projectName == projInfo.projName:
                return True
        return False
    
        
        