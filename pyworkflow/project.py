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
This modules handles the Project management
"""

import os
from os.path import abspath, split

from pyworkflow.em import *
from pyworkflow.apps.config import *
from pyworkflow.protocol import *
from pyworkflow.mapper import SqliteMapper
from pyworkflow.utils import cleanPath, makePath, makeFilePath, join, exists, runJob, copyFile
from pyworkflow.hosts import HostMapper, HostConfig
import pyworkflow.protocol.launch as jobs

PROJECT_DBNAME = 'project.sqlite'
PROJECT_LOGS = 'Logs'
PROJECT_RUNS = 'Runs'
PROJECT_TMP = 'Tmp'
PROJECT_SETTINGS = 'settings.sqlite'

class Project(object):
    """This class will handle all information 
    related with a Project"""
    def __init__(self, path):
        """For create a Project, the path is required""" 
        self.name = path
        self.path = abspath(path)
        self.pathList = [] # Store all related paths
        self.dbPath = self.addPath(PROJECT_DBNAME)
        self.logsPath = self.addPath(PROJECT_LOGS)
        self.runsPath = self.addPath(PROJECT_RUNS)
        self.tmpPath = self.addPath(PROJECT_TMP)
        self.settingsPath = self.addPath(PROJECT_SETTINGS)
        self.runs = None
        
    def getObjId(self):
        """ Return the unique id assigned to this project. """
        return os.path.basename(self.path)
    
    def addPath(self, *paths):
        """Store a path needed for the project"""
        p = self.getPath(*paths)
        self.pathList.append(p)
        return p
        
    def getPath(self, *paths):
        """Return path from the project root"""
        return join(*paths)
    
    def getName(self):
        return self.name
    
    def getTmpPath(self, *paths):
        return self.getPath(PROJECT_TMP, *paths)
    
    def getSettings(self):
        return self.settings
    
    def saveSettings(self):
        self.settings.write()
    
    def load(self):
        """Load project data and settings
        from the project dir."""
        
        if not exists(self.path):
            raise Exception("Project doesn't exists in '%s'" % self.path)
        os.chdir(self.path) #Before doing nothing go to project dir
        if not exists(self.dbPath):
            raise Exception("Project database not found in '%s'" % self.dbPath)
        self.mapper = SqliteMapper(self.dbPath, globals())
        self.settings = loadSettings(self.settingsPath)
        #self.hostsMapper = HostMapper(self.settingsPath)
        
    def create(self, defaultSettings):
        """Prepare all required paths and files to create a new project.
        Params:
         hosts: a list of configuration hosts associated to this projects (class ExecutionHostConfig)
        """
        #cleanPath(self.settingsPath)
        # Create project path if not exists
        makePath(self.path)
        os.chdir(self.path) #Before doing nothing go to project dir
        self.clean()
        print abspath(self.dbPath)
        # Create db throught the mapper
        self.mapper = SqliteMapper(self.dbPath, globals())
        self.mapper.commit()
        # Write settings to disk
        self.settings = defaultSettings
        self.settings.write(self.settingsPath)
        
        # Create other paths inside project
        for p in self.pathList:
            if '.' in p:
                makeFilePath(p)
            else:
                makePath(p)
        
    def clean(self):
        """Clean all project data"""
        cleanPath(*self.pathList)      
                
    def launchProtocol(self, protocol, wait=False):
        """ In this function the action of launching a protocol
        will be initiated. Actions done here are:
        1. Store the protocol and assign name and working dir
        2. Create the working dir and also the protocol independent db
        3. Call the launch method in protocol.job to handle submition: mpi, thread, queue,
        and also take care if the execution is remotely."""
        #print ">>> PROJECT: launchProtocol"
        protocol.setStatus(STATUS_LAUNCHED)
        self._setupProtocol(protocol)
        
        #protocol.setMapper(self.mapper) # mapper is used in makePathAndClean
        protocol.makePathsAndClean() # Create working dir if necessary
        self.mapper.commit()
        
        # Prepare a separate db for this run
        # NOTE: now we are simply copying the entire project db, this can be changed later
        # to only create a subset of the db need for the run
        copyFile(self.dbPath, protocol.getDbPath())
        
        # Launch the protocol, the jobId should be set after this call
        jobs.launch(protocol, wait)
        
        # Commit changes
        if wait: # This is only useful for launching tests...
            self._updateProtocol(protocol)
        else:
            self.mapper.store(protocol)
        self.mapper.commit()
        #print ">>> PROJECT: launchProtocol: DONE"
        
    def _updateProtocol(self, protocol):
        try:
            # FIXME: this will not work for a real remote host
            jobId = protocol.getJobId() # Preserve the jobId before copy
            
            dbPath = self.getPath(protocol.getDbPath())
            #join(protocol.getHostConfig().getHostPath(), protocol.getDbPath())
            prot2 = getProtocolFromDb(dbPath, protocol.getObjId(), globals())
            # Copy is only working for db restored objects
            protocol.setMapper(self.mapper)
            protocol.copy(prot2)
            #protocol.mapper.copyRelations(prot2.mapper.getRelations(prot2))
            # Restore jobId
            protocol.setJobId(jobId)
            
            self.mapper.store(protocol)
        except Exception, ex:
            print "Error trying to update protocol: %s\n %s" % (jobId, ex)
            
        
    def stopProtocol(self, protocol):
        """ Stop a running protocol """
        jobs.stop(protocol)
        protocol.setAborted()
        self._storeProtocol(protocol)
        
    def continueProtocol(self, protocol):
        """ This function should be called 
        to mark a protocol that have an interactive step
        waiting for approval that can continue
        """
        for step in protocol._steps:
            if step.status == STATUS_WAITING_APPROVAL:
                step.setStatus(STATUS_FINISHED)
                self.mapper.store(step)
                self.mapper.commit()
                break
        self.launchProtocol(protocol)
        
    def deleteProtocol(self, protocol):
        self.mapper.delete(protocol) # Delete from database
        wd = protocol.workingDir.get()
        if wd.startswith(PROJECT_RUNS):
            cleanPath()
        else:
            print "Error path: ", wd 
      
        self.mapper.commit()     
        
    def copyProtocol(self, protocol):
        """ Make a copy of the protocol, return a new one with copied values. """
        cls = protocol.getClass()
        newProt = cls() # Create new protocol instance
        newProt.copyDefinitionAttributes(protocol)
        newProt.setMapper(self.mapper)
        
        return newProt
    
    def saveProtocol(self, protocol):
        protocol.setStatus(STATUS_SAVED)
        if protocol.hasObjId():
            self._storeProtocol(protocol)
        else:
            self._setupProtocol(protocol)
        
    def _setHostConfig(self, protocol):
        """ Set the appropiate host config to the protocol
        give its value of 'hostname'
        """
        hostName = protocol.getHostName()
        hostConfig = self.settings.getHostByLabel(hostName)
        hostConfig.cleanObjId()
        # Add the project name to the hostPath in remote execution host
        hostRoot = hostConfig.getHostPath()
        hostConfig.setHostPath(join(hostRoot, os.path.basename(self.path)))
        protocol.setHostConfig(hostConfig)
    
    def _storeProtocol(self, protocol):
        self.mapper.store(protocol)
        self.mapper.commit()
                
    def _setupProtocol(self, protocol):
        """Insert a new protocol instance in the database"""
        self._storeProtocol(protocol) # Store first to get a proper id
        # Set important properties of the protocol
        name = protocol.getClassName() + protocol.strId()
        protocol.setName(name)
        protocol.setWorkingDir(self.getPath(PROJECT_RUNS, name))
        protocol.setMapper(self.mapper)
        self._setHostConfig(protocol)
        # Update with changes
        self._storeProtocol(protocol)
        
    def getRuns(self, iterate=False, refresh=True):
        """ Return the existing protocol runs in the project. 
        """
        if self.runs is None or refresh:
            self.runs = self.mapper.selectByClass("Protocol", iterate=False)
            for r in self.runs:
                if r.isActive():
                    self._updateProtocol(r)
            self.mapper.commit()
        
        return self.runs
    
    def getRunsGraph(self, refresh=True):
        """ Build a graph taking into account the dependencies between
        different runs, ie. which outputs serves as inputs of other protocols. 
        """
        #import datetime as dt # TIME PROFILE
        #t = dt.datetime.now()
        outputDict = {} # Store the output dict
        runs = self.getRuns(refresh=refresh)
        from pyworkflow.utils.graph import Graph
        g = Graph(rootName='PROJECT')
        
        for r in runs:
            n = g.createNode(r.strId())
            n.run = r
            n.label = r.getRunName()
            
            for _, attr in r.iterOutputAttributes(EMObject):
                outputDict[attr.getName()] = n # mark this output as produced by r
            
        def _checkInputAttr(node, pointed):
            """ Check if an attr is registered as output"""
            pointedName = pointed.getName()
            if pointedName in outputDict:
                parentNode = outputDict[pointedName]
                parentNode.addChild(node)
                return True
            return False
            
        for r in runs:
            node = g.getNode(r.strId())
            for _, attr in r.iterInputAttributes():
                if attr.hasValue():
                    pointed = attr.get()
                    # Only checking pointed object and its parent, if more levels
                    # we need to go up to get the correct dependencies
                    _checkInputAttr(node, pointed) or _checkInputAttr(node, self.mapper.getParent(pointed))
                    
        rootNode = g.getRoot()
        rootNode.run = None
        rootNode.label = "PROJECT"
        
        for n in g.getNodes():
            if n.isRoot() and not n is rootNode:
                rootNode.addChild(n)
        return g
        
