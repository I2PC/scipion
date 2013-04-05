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
from pyworkflow.mapper import SqliteMapper
from pyworkflow.utils.path import cleanPath, makePath, getHomePath


PROJECTS_PATH = 'Scipion_Projects'


class Manager(object):
    """This class will handle the creation, modification
    and listing of projects."""
    def __init__(self):
        """For create a Project, the path is required"""
        self.path = join(getHomePath(), PROJECTS_PATH)
        makePath(self.path)
        
    def getProjectPath(self, projectName):
        """Return the project path given the name"""
        return join(self.path, projectName)
        
    def listProjects(self):
        """Return a list with all existing projects"""
        return [f for f in os.listdir(self.path) if os.path.isdir(self.getProjectPath(f))]
    
    def createProject(self, projectName):
        """Create a new project """
        proj = Project(self.getProjectPath(projectName))
        proj.create()
        
    def deleteProject(self, projectName):
        pass
        
        