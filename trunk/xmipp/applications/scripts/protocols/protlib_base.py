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
    def __init__(self, projectDir):
        self.projectDir = projectDir
        self.cfgName = projectDefaults['Cfg']
        self.dbName =  projectDefaults['Db']
        self.logsDir = projectDefaults['LogsDir']
        self.runsDir = projectDefaults['RunsDir']
    
    def exists(self):
        ''' a project exists if the data base can be opened and directories Logs and Runs
        exists'''
        from protlib_sql import existsDB
        status = True
        if not os.path.exists(self.logsDir):
            status = False
        elif not os.path.exists(self.runsDir):
            status = False
        elif not os.path.exists(self.cfgName):
            status = False
        elif not existsDB(self.dbName):
            status = False
        return status
    
    def create(self):
        os.chdir(self.projectDir)
        print "Creating new project on directory: '%s'" % self.projectDir
         #==== CREATE CONFIG file
        self.config = ConfigParser.RawConfigParser()            
        self.config.add_section('project')
        self.config.set('project', 'projectdir', self.projectDir)
        self.config.set('project', 'mpiflavour', 'OPEN-MPI')
        self.writeConfig()
        #===== CREATE LOG AND RUN directories
        if not os.path.exists(self.logsDir):
            os.mkdir(self.logsDir)
        if not os.path.exists(self.runsDir):
            os.mkdir(self.runsDir)
        #===== CREATE DATABASE
        self.projectDb  = XmippProjectDb(self.dbName)
        #===== POPULATE SOME TABLES
        for section, groupList in sections:
            for group in groupList:
                groupName = group[0]            
                prots = group[1:]
                self.projectDb.insertGroup(groupName)
                for p in prots:
                    self.projectDb.insertProtocol(groupName, launchDict[p])
        # commint changes
        self.projectDb.connection.commit()
        
    def load(self):
        print 'Loading project..'
        self.config = ConfigParser.RawConfigParser()
        self.config.read(self.cfgName)
        # Load database
        self.projectDb = XmippProjectDb(self.dbName)
        
    def writeConfig(self):
       with open(self.cfgName, 'wb') as configfile:
            self.config.write(configfile) 
            

class XmippProtocol(object):
    '''This class will serve as base for all Xmipp Protocols'''
    
    def __init__(self, protocolName, runName, project=None):
        '''Basic constructor of the Protocol
        protocolName -- the name of the protocol, should be unique
        runName      -- the name of the run,  should be unique for one protocol
        project      -- project instance
        '''

        self.Name = protocolName
        self.runName = runName
        #A protocol must be able to find its own project
        self.project = project
        self.Import = '' # this can be used by database for import modules
        self.WorkingDir = os.path.join(launchDict['Projection Matching'],runName)
        self.ProjectDir = project.ProjectDir  
        #Setup the Log for the Protocol
        self.LogDir = project.logsDir
        uniquePrefix = self.WorkingDir.replace('/', '_')
        self.LogPrefix = os.path.join(self.LogDir, uniquePrefix)       
        self.errors = []
        self.summary = []
        self.continueAt=1
        self.isIter=False
        
    def getProjectId(self):
        pass
        #self.project = project.getId(launchDict['Projection Matching'],runName,)
        
    def validate(self):
        '''Validate if the protocols is ready to be run
        it may be redefined in derived protocol classes but do not forget
        to call the main class with
        super(ProtProjMatch, self).validate()
        '''
        self.errors=[]
        #check if there is a valid project, otherwise abort
        if not self.project.exists():
            self.errors.append("Not Valid project available")
    
    def summary(self):
        '''Produces a summary with the most relevant information of the protocol run'''
        self.summary=[]
        self.summary.append(self.Name)
    
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
        #self.project.load()
        self.Db = XmippProtocolDb(project.dbName, self.continueAt,self.isIter,self.run_id)

    def __init__(self, continueAt, isIter, run_id):

        self.init()
        
    def command_line_options(self):
        '''process protocol command line'''
        import optparse
        self.parser = optparse.OptionParser()
        self.parser.add_option('-g', '--gui',
                                  dest="gui",
                                  default=False,
                                  action="store_true",
                                  help="use graphic interface to launch protocol "
                                  )
        self.parser.add_option('-c', '--no_check',
                                  dest="no_check",
                                  default=False,
                                  action="store_true",
                                  help="do NOT check run checks before execute protocols"
                                  )
        self.gui, self.no_check = self.parser.parse_args()

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
    
    
      

    
    

            
      

