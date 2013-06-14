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
from pyworkflow.utils import cleanPath, makePath, makeFilePath, join, exists, runJob
from pyworkflow.apps.config import ExecutionHostMapper, ExecutionHostConfig
from pyworkflow.protocol.launch import launchProtocol

PROJECT_DBNAME = 'project.sqlite'
PROJECT_LOGS = 'Logs'
PROJECT_RUNS = 'Runs'
PROJECT_HOSTS = 'execution_hosts.xml'

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
        self.hostsPath = self.addPath(PROJECT_HOSTS)
        
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
    
    def load(self):
        """Load project data and settings
        from the project dir."""
        
        if not exists(self.path):
            raise Exception("Project doesn't exists in '%s'" % self.path)
        os.chdir(self.path) #Before doing nothing go to project dir
        if not exists(self.dbPath):
            raise Exception("Project database not found in '%s'" % self.dbPath)
        self.mapper = SqliteMapper(self.dbPath, globals())
        self.hostsMapper = ExecutionHostMapper(self.hostsPath)
        
    def create(self, hosts):
        """Prepare all required paths and files to create a new project.
        Params:
         hosts: a list of configuration hosts associated to this projects (class ExecutionHostConfig)
        """
        #cleanPath(self.hostsPath)
        # Create project path if not exists
        makePath(self.path)
        os.chdir(self.path) #Before doing nothing go to project dir
        self.clean()
        print abspath(self.dbPath)
        # Create db throught the mapper
        self.mapper = SqliteMapper(self.dbPath, globals())
        self.mapper.commit()
        # Write hosts configuration to disk
        self.hostsMapper = ExecutionHostMapper(self.hostsPath)

        for h in hosts:
            self.hostsMapper.insert(h)
        self.hostsMapper.commit()
        # Create other paths inside project
        for p in self.pathList:
            print "creating file: ", p
            if '.' in p:
                makeFilePath(p)
            else:
                makePath(p)
        
    def clean(self):
        """Clean all project data"""
        cleanPath(*self.pathList)      
                
    def launchProtocol(self, protocol, wait=False):
        """Launch another scritp to run the protocol."""
        # TODO: Create a launcher class that will 
        # handle the communication of remote projects
        # and also the particularities of job submission: mpi, threads, queue, bash
        self._setupProtocol(protocol)
        jobId = launchProtocol(protocol, wait)
        if wait:
            self.mapper.updateFrom(protocol)
        
        
    def runProtocol(self, protocol, mpiComm=None):
        """Directly execute the protocol.
        mpiComm is only used when the protocol steps are execute
        with several MPI nodes.
        """
        protocol.mapper = self.mapper
        protocol._stepsExecutor = StepExecutor()
        if protocol.stepsExecutionMode == STEPS_PARALLEL:
            if protocol.numberOfMpi > 1:
                protocol._stepsExecutor = MPIStepExecutor(protocol.numberOfMpi.get(), mpiComm)
            elif protocol.numberOfThreads > 1:
                protocol._stepsExecutor = ThreadStepExecutor(protocol.numberOfThreads.get()) 
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
        
    def copyProtocol(self, protocol):
        """ Make a copy of the protocol, return a new one with copied values. """
        cls = protocol.getClass()
        newProt = cls() # Create new protocol instance
        newProt.copyDefinitionAttributes(protocol)
        
        return newProt
    
    def saveProtocol(self, protocol):
        protocol.status.set(STATUS_SAVED)
        self._storeProtocol(protocol)
        
    def _setHostConfig(self, protocol):
        """ Set the appropiate host config to the protocol
        give its value of 'hostname'
        """
        hostName = protocol.getHostName()
        print "hostName: ", hostName
        hostConfig = self.hostsMapper.selectByLabel(hostName)
        hostConfig.cleanObjId()
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
        protocol.setMapper(self.mapper)
        protocol.dbPath = self.dbPath
        protocol.setWorkingDir(self.getPath("Runs", name))
        self._setHostConfig(protocol)
        # Update with changes
        self._storeProtocol(protocol)
        
    def getHosts(self):
        """ Retrieve the hosts associated with the project. (class ExecutionHostConfig) """
        return self.hostsMapper.selectAll()
    
    def launchRemoteProtocol(self, protocol):
        """ Launch protocol in an execution host    
        Params:
            potocol: Protocol to launch in an execution host.
        """
        # Fisrt we must recover the execution host credentials.
        self.hostsMapper = ExecutionHostMapper(self.hostsPath)
        executionHostConfig = self.hostsMapper.selectByLabel(protocol.getHostName())
        self.sendProtocol(protocol, executionHostConfig)
        command = "nohup pw_protocol_run.py " + self.getObjId() + " " + str(protocol.getObjId()) + " > /gpfs/fs1/home/bioinfo/apoza/salida.txt 2> /gpfs/fs1/home/bioinfo/apoza/salida2.txt &"
        from pyworkflow.utils.utils import executeRemote
        executeRemote(command,
                      executionHostConfig.getHostName(),
                      executionHostConfig.getUserName(),
                      executionHostConfig.getPassword())
        
    def sendProtocol(self, protocol, executionHostConfig):
        """ Send protocol to an execution host    
        Params:
            potocol: Protocol to send to an execution host.
        """
        from utils.file_transfer import FileTransfer
        filePathDict = {}
        # We are going to create project folder in the remote host
        projectFolder = split(self.path)[1]
        # Prepare source and target files
        for filePath in protocol.getFiles():
            sourceFilePath = os.path.join(self.path, filePath)
            targetFilePath = join(executionHostConfig.getHostPath(), projectFolder, filePath)
            filePathDict[sourceFilePath] = targetFilePath
            
        # We add sqlite database file
        filePathDict[join(self.path, self.dbPath)] = join(executionHostConfig.getHostPath(), projectFolder, self.dbPath)
        # Transfer files       
        fileTransfer = FileTransfer()
        fileTransfer.transferFilesTo(filePathDict,
                        executionHostConfig.getHostName(),
                        executionHostConfig.getUserName(),
                        executionHostConfig.getPassword(),
                        gatewayHosts = None, 
                        numberTrials = 1,                        
                        forceOperation = False,
                        operationId = 1)
        
