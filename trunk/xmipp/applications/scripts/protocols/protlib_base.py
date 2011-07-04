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
import sys

class XmippProject():
    def __init__(self, projectDir=None):
        if not projectDir:
            self.projectDir = os.getcwd()
        else:
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
        # commit changes
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
            
    def clean(self):
        ''' delete files related with a project'''
        if not os.path.exists(self.logsDir):
            shutil.rmtree(self.logsDir)
        if not os.path.exists(self.runsDir):
            shutil.rmtree(self.runsDir)
        if not os.path.exists(self.cfgName):
            os.remove(self.cfgName)
        if not os.path.exists(self.dbName):
            os.remove(self.dbName)
            
class XmippProtocol(object):
    '''This class will serve as base for all Xmipp Protocols'''
    def __init__(self, protocolName, scriptname,runName, project=None, NumberOfMpiProcesses=None):
        '''Basic constructor of the Protocol
        protocolName -- the name of the protocol, should be unique
        scriptname,     file containing the protocol instance
        runName      -- the name of the run,  should be unique for one protocol
        project      -- project instance
        '''
        self.DoParallel = True
        try:
            NumberOfMpiProcesses
            if NumberOfMpiProcesses == None or NumberOfMpiProcesses ==0:
                self.DoParallel = False
        except NameError:
            self.DoParallel = False
        #do something similar with mpiflavour whenever we decide what to do
#        self.mpiflavour = ''
#        try:
#            SystemFlavour
#            if SystemFlavour == None or SystemFlavour !='':
#                self.mpiflavour = SystemFlavour
#        except NameError:
#            self.DoParallel = False
#        SystemFlavour
            
        self.Name = protocolName
        self.runName = runName
        self.scriptName = scriptname
        #A protocol must be able to find its own project
        self.project = project
        self.Import = '' # this can be used by database for import modules
        self.WorkingDir = os.path.join(protocolName,runName)
        self.projectDir = project.projectDir  
        #Setup the Log for the Protocol
        self.LogDir = project.logsDir
        self.uniquePrefix = self.WorkingDir.replace('/', '_')
        self.LogPrefix = os.path.join(self.LogDir, self.uniquePrefix)       
        self.errors = []
        self.summary = []
        self.continueAt=1
        self.isIter=False
        self.DoDeleteWorkingDir=False
        
        
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
            self.errors.append("Not valid project available")
    
    def summary(self):
        '''Produces a summary with the most relevant information of the protocol run'''
        self.summary=[]
        self.summary.append(self.comment)
    
    def warnings(self):
        '''Output some warnings that can be errors and require user confirmation to procceed'''
        pass
    
    def postRun(self):
        '''This function will be called after the run function is executed'''
        pass   
    
    def defineActions(self):
        '''In this function the actions to be performed by the protocol will be added to db.
        each particular protocol need to add its specific actions. Thing to be performed before 
        "run" should be added here'''
        pass
    
    def runSetup(self):
        logfile  = self.LogPrefix + ".log"
        self.Log = XmippLog(logfile, logfile)
        self.Db  = XmippProtocolDb(self.project.dbName
                                , self.restartStep
                                , self.isIter
                                , self)

        

    def run(self, restartStep=1, isIter=True):
        '''Run of the protocols
        if the other functions have been correctly implemented, this not need to be
        touched in derived class, since the run of protocols should be the same'''
        
        #Change to project dir
        os.chdir(self.projectDir);

        errors = self.validate()
        if len(errors) > 0:
            raise Exception('\n'.join(errors))
        self.restartStep = restartStep
        self.isIter = isIter
        #Initialization of log and db
        self.runSetup()
        
        self.Db.setIteration(1)
        #insert basic operations for all scripta
        self.Db.insertAction('deleteWorkingDirectory'
                                                        ,DoDeleteWorkingDir  = self.DoDeleteWorkingDir
                                                        ,WorkingDir          = self.WorkingDir)

        self.Db.insertAction('createDir', path = self.WorkingDir)
        self.Db.setIteration(XmippProtocolDbStruct.doAlways)
        self.Db.insertAction('makeScriptBackup',script     = self.scriptName
                                                      ,WorkingDir = self.WorkingDir)#backup always
        self.Db.setIteration(1)

        #Add actions to database
        self.defineActions()
        #Run actions from database
        self.Db.runActions(self.Log, self.Import)
        #Stuff after running
        self.postRun()
        return 0
    
def command_line_options():
    '''process protocol command line'''
    try:
        import argparse
        parser = argparse.ArgumentParser(description='test goes here')
        parser.add_argument('-g', '--gui',
                                  dest="gui",
                                  default=False,
                                  action="store_true",
                                  help="use graphic interface to launch protocol "
                                  )
        parser.add_argument('-c', '--no_check',
                                  dest="no_check",
                                  default=False,
                                  action="store_true",
                                  help="do NOT check run checks before execute protocols"
                                  )
        (options, args) = parser.parse_args()
        
    except:
        import optparse
        parser = optparse.OptionParser()
        parser.add_option('-g', '--gui',
                                  dest="gui",
                                  default=False,
                                  action="store_true",
                                  help="use graphic interface to launch protocol "
                                  )
        parser.add_option('-c', '--no_check',
                                  dest="no_check",
                                  default=False,
                                  action="store_true",
                                  help="do NOT check run checks before execute protocols"
                                  )
        (options, args) = parser.parse_args()
    return options

def protocolMain(ProtocolClass):
    script  = sys.argv[0]
    options = command_line_options()
    
    #1) just call the gui
    
    if options.gui:
        from protlib_gui import ProtocolGUI 
        gui = ProtocolGUI()
        gui.createGUI(script)
        gui.fillGUI()
        gui.launchGUI()
    else:#2) Run from command line
        mod = loadModule(script)
        #init project
        project = XmippProject()
        #load project: read config file and open conection database
        project.load()
        #register run with runName, script, comment=''):
        p = ProtocolClass(script, project)
        run_id = project.projectDb.getRunId(p.Name, mod.RunName)
        
        if options.no_check: 
            if not run_id:
                reportError("Protocol run '%s' has not been registered in project database" % mod.RunName)
        else:
            if not confirmWarning(p.warnings()):
                exit(0)
            
            if not run_id:
                _run={
                       'comment' : mod.Comment
                      ,'protocol_name': p.Name
                      ,'run_name' : mod.RunName
                      ,'script'  : script 
                      }
                project.projectDb.insertRun(_run)
            if mod.SubmmitToQueue:
                submitProtocol(script,
                               jobId = p.uniquePrefix,
                               queueName = mod.QueueName,
                               nodes = mod.NumberOfMpiProcesses,
                               threads = mod.NumberOfThreads,
                               hours = mod.QueueHours,
                               command = 'python %s --no_check' % script
                               )
                exit(0)
#            else:
#                reportError("%s[ Protocol runs already created, data will be overwrited ]%s\n"  << quitar color 
#                                 % (bcolors.FAIL, bcolors.ENDC)) 
        
        p.run(mod.ContinueAtIteration,mod.IsIter)


