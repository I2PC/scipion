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
import log
import logging
from protlib_sql import *
from protlib_utils import *

def checkErrors():
    '''This function will be used to validate the protocols
    should be implemented in all derived protocols'''
    return []

class XmippProtocol(object):
    '''This class will serve as base for all Xmipp Protocols'''
    
    def __init__(self, scriptname, workingdir, projectdir=None, logdir='Logs', restartStep=1, isIter=True):
        '''Basic constructor of the Protocol'''
        '''
        scriptname  -- the name of the protocol script name, ej: xmipp_protocol_ml2d.py
        projectdir  -- directory of the project, usually where to run several protocols
        workingdir  -- directory for the output of this protocol (relative to projectdir) 
        
        logdir      -- directory for logs (relative to projectdir)
        restartStep -- at wich step do you wish to continue the protocol if was launched previously
        isIter      -- if True previous param refers to iter number, if not, to step number 
        '''
        self.Name = scriptname
        #Setup prefix from scriptname, this impose using names 'xmipp_protocol_xxx'
        self.Prefix = getScriptPrefix(scriptname)
        self.Import = '' # this can be used by database for import modules
        self.WorkingDir = workingdir        
        self.Step = restartStep        
                
        if not projectdir:
            self.ProjectDir = os.getcwd()
        else:
            self.ProjectDir = projectdir
                
        #Setup the Log for the Protocol
        self.LogDir = logdir
        uniquePrefix = "%s_%s" % (self.Prefix, workingdir.replace('/', '_'))
        self.LogPrefix = os.path.join(logdir, uniquePrefix)

        #Create dir if not exists
        if not os.path.exists(logdir):
            os.makedirs(logdir)            
        
        logfile = self.LogPrefix + ".log"
        self.Log = XmippLog( logfile, logfile )
        #Setup database for executing commands
        dbfile = self.LogPrefix + ".sqlite"
        self.Db = XmippProtocolDb(dbfile, self.Prefix + "Table", self.Step, isIter)

        
    def validate(self):
        '''This function will validate if the protocols is ready to be run
        it should be redefine in derived protocol classes, it will be a wrapper
        around the module function checkErrors
        '''
        return checkErrors()
    
    def preRun(self):
        '''This function will be called before the run function is executed'''
        pass
    
    def defineActions(self):
        '''In this funciton the actions to be performed by the protocol will be add to db.
        each particular protocol need to add its specific actions'''
        pass
    
    def run(self):
        '''Run of the protocols
        if the other functions have been correctly implemented, this not need to be
        touched in derived class, since the run of protocols should be the same'''
        errors = self.validate()
        if len(errors) > 0:
            raise Exception('\n'.join(errors))
        #Stuff before running
        self.preRun()
        #Add actions to database
        self.defineActions()
        #Change to project dir
        os.chdir(self.ProjectDir);
        #Run actions from database
        self.Db.runActions(self.Log, self.Import)
        return 0
      
class XmippProject():
    pass
            
      

