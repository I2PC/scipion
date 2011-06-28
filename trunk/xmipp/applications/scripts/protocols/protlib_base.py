#!/usr/bin/env python
'''
#/***************************************************************************
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
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
# ***************************************************************************
 '''
 
import os
import shutil
import ConfigParser
from config import *
from protlib_sql import *
from protlib_utils import *

class XmippProject():
    def __init__(self, dir):
        self.dir = dir
        self.cfgName = projectDefaults['Cfg']
        self.dbName =  projectDefaults['Db']
        self.logsDir = projectDefaults['LogsDir']
        self.runsDir = projectDefaults['RunsDir']
    
    def exists(self):
        pass
    
    def create(self):
        os.chdir(self.dir)
        print "Creating new project on directory: '%s'" % self.dir
         #==== CREATE CONFIG file
        self.config = ConfigParser.RawConfigParser()            
        self.config.add_section('project')
        self.config.set('project', 'projectdir', self.dir)
        self.config.set('project', 'mpiflavour', 'OPEN-MPI')
        self.writeConfig()
        #===== CREATE LOG AND RUN directories
        os.mkdir(self.logsDir)
        os.mkdir(self.runsDir)
        #===== CREATE DATABASE
        self.db = XmippProjectDb(self.dbName)
        #===== POPULATE SOME TABLES
        for section, groupList in sections:
            for group in groupList:
                groupName = group[0]            
                prots = group[1:]
                self.db.insertGroup(groupName)
                for p in prots:
                    self.db.insertProtocol(groupName, launchDict[p])
        # commint changes
        self.db.connection.commit()
        
    def load(self):
        print 'Loading project..'
        self.config = ConfigParser.RawConfigParser()
        self.config.read(self.cfgName)
        # Load database
        self.db = XmippProjectDb(self.dbName)
        
    def writeConfig(self):
       with open(self.cfgName, 'wb') as configfile:
            self.config.write(configfile) 
            

class XmippProtocol(object):
    '''This class will serve as base for all Xmipp Protocols'''
    
    def __init__(self, protocolName, runName, project):
        '''Basic constructor of the Protocol
        protocolName -- the name of the protocol, should be unique
        runName      -- the name of the run,  should be unique for one protocol
        project      -- project instance
        '''
        self.Name = protocolName
        self.runName = runName
        self.project = project
        self.Import = '' # this can be used by database for import modules
        self.WorkingDir = os.path.join(protocolName, runName)
        self.ProjectDir = project.dir  
        #Setup the Log for the Protocol
        self.LogDir = project.logsDir
        uniquePrefix = self.WorkingDir.replace('/', '_')
        self.LogPrefix = os.path.join(self.LogDir, uniquePrefix)       

        
    def validate(self):
        '''Validate if the protocols is ready to be run
        it should be redefine in derived protocol classes, it will be a wrapper
        around the module function checkErrors
        '''
        pass
    
    def summary(self):
        '''Produces a summary with the most relevant information of the protocol run'''
        pass
    
    def warnings(self):
        '''Output some warnings that can be errors and require user confirmation to procceed'''
        pass
    
    def preRun(self):
        '''This function will be called before the run function is executed'''
        pass
    
    def postRun(self):
        '''This function will be called before the run function is executed'''
        pass   
    
    def defineActions(self):
        '''In this funciton the actions to be performed by the protocol will be add to db.
        each particular protocol need to add its specific actions'''
        pass
    
    def init(self):
        #Create dir if not exists
        #if not os.path.exists(logdir):
        #    os.makedirs(logdir)            
        
        logfile = self.LogPrefix + ".log"
        self.Log = XmippLog(logfile, logfile)
        #Setup database for executing commands
        #dbfile = self.LogPrefix + ".sqlite"
        #self.Db = self.project.XmippProtocolDb(dbfile, self.Prefix + "Table", self.Step, isIter)
        self.Db = project.db
        
    def run(self, restartStep=1, isIter=True):
        '''Run of the protocols
        if the other functions have been correctly implemented, this not need to be
        touched in derived class, since the run of protocols should be the same'''
        
        errors = self.validate()
        if len(errors) > 0:
            raise Exception('\n'.join(errors))
        self.Step = restartStep
        #Initialization of log and db
        self.init()        
        #Add actions to database
        self.defineActions()
        #Change to project dir
        os.chdir(self.ProjectDir);
        #Stuff before running
        self.preRun()
        #Run actions from database
        self.Db.runActions(self.Log, self.Import)
        #Stuff after running
        self.postRun()
        return 0
    
    
      

    
    

            
      

