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
from protocol.protocol import STATUS_WAITING_APPROVAL, STATUS_FINISHED
"""
This modules handles the Project management
"""
import os
from os.path import abspath

from pyworkflow.mapper import SqliteMapper
from pyworkflow.utils.path import cleanPath, makePath, join, exists
from pyworkflow.protocol import *
from pyworkflow.em import *


PROJECT_DBNAME = 'project.sqlite'
PROJECT_LOGS = 'Logs'
PROJECT_RUNS = 'Runs'


class Project(object):
    """This class will handle all information 
    related with a Project"""
    def __init__(self, path):
        """For create a Project, the path is required""" 
        self.path = abspath(path)
        self.pathList = [] # Store all related paths
        self.dbPath = self.addPath(PROJECT_DBNAME)
        self.logsPath = self.addPath(PROJECT_LOGS)
        self.runsPath = self.addPath(PROJECT_RUNS)    
    
    def addPath(self, *paths):
        """Store a path needed for the project"""
        p = self.getPath(*paths)
        self.pathList.append(p)
        return p
        
    def getPath(self, *paths):
        """Return path from the project root"""
        return join(*paths)
    
    def load(self):
        """Load project data and settings
        from the project dir."""
        os.chdir(self.path) #Before doing nothing go to project dir
        if not exists(self.dbPath):
            raise Exception("Project doesn't exists in '%s'" % self.path)
        self.mapper = SqliteMapper(self.dbPath, globals())
        
    def create(self):
        """Create a new project.
        Prepare all required paths and files"""
        # Create project path if not exists
        makePath(self.path)
        os.chdir(self.path) #Before doing nothing go to project dir
        print abspath(self.dbPath)
        # Create db throught the mapper
        self.mapper = SqliteMapper(self.dbPath, globals())
        self.mapper.commit()
        # Create other paths inside project
        makePath(*self.pathList)
        
    def clean(self):
        """Clean all project data"""
        cleanPath(*self.pathList)        
                
    def launchProtocol(self, protocol, wait=False):
        """Launch another scritp to run the protocol."""
        # TODO: Create a launcher class that will 
        # handle the communication of remove projects
        # and also the particularities of job submission: mpi, threads, queue, bash
        self.insertProtocol(protocol)        
        from pyworkflow.protocol.launcher import FULLPATH
        cmd = "python %s %s %s" % (FULLPATH, os.path.basename(self.path), 
                                     protocol.strId())
        if not wait:
            cmd += ' &'
        print cmd
        os.system(cmd)
        
    def runProtocol(self, protocol):
        """Directly execute the protocol"""
        protocol.mapper = self.mapper
        protocol.run()
        
    def continueProtocol(self, protocol):
        """ This function should be called 
        to mark a protocol that have an interactive step
        waiting for approval that can continue
        """
        for step in protocol._steps:
            if step.status == STATUS_WAITING_APPROVAL:
                step.status.set(STATUS_FINISHED)
                self.mapper.store(step)
                self.mapper.commit()
                break
        self.launchProtocol(protocol)
        
    def deleteProtocol(self, protocol):
        self.mapper.delete(protocol) # Delete from database
        cleanPath(protocol.workingDir.get())  
        self.mapper.commit()     
        
    def insertProtocol(self, protocol):
        """Insert a new protocol instance in the database"""
        self.mapper.store(protocol)
        name = protocol.getClassName() + protocol.strId()
        protocol.setName(name)
        protocol.mapper =  self.mapper
        protocol.workingDir.set(self.getPath("Runs", name))
        protocol._store() #FIXME, actually only need to update workingDir and name
        
        