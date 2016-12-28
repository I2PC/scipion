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
This modules contains classes to store information about
execution hosts.
"""

import sys
import json
from ConfigParser import ConfigParser
from collections import OrderedDict

import pyworkflow as pw
from pyworkflow.object import *
from pyworkflow.mapper import SqliteMapper, XmlMapper



class HostMapper(SqliteMapper):
    def __init__(self, filename, dictClasses=None):
        if dictClasses is None:
            dictClasses = globals()
        SqliteMapper.__init__(self, filename, dictClasses)
        
    def selectByLabel(self, objLabel):
        hostsList = self.selectAll()
        for host in hostsList:
            if host.label == objLabel:
                return host
        return None
        
        
class HostConfig(OrderedObject):
    """ Main store the configuration for execution hosts. """
    
    def __init__(self, **kwargs):
        OrderedObject.__init__(self, **kwargs)
        self.label = String(kwargs.get('label', None))
        self.hostName = String(kwargs.get('hostName', None))
        self.userName = String()
        self.password = String()
        self.hostPath = String()
        self.mpiCommand = String()
        self.scipionHome = String()
        self.scipionConfig = String()
        self.address = String()
        self.queueSystem = QueueSystemConfig()
    
    def getLabel(self):
        return self.label.get()
    
    def getHostName(self):
        return self.hostName.get()
    
    def getUserName(self):
        return self.userName.get()
    
    def getPassword(self):
        return self.password.get()
    
    def getHostPath(self):
        return self.hostPath.get()
    
    def getSubmitCommand(self):
        return self.queueSystem.submitCommand.get()
    
    def getSubmitPrefix(self):
        return self.queueSystem.submitPrefix.get()
    
    def getCancelCommand(self):
        return self.queueSystem.cancelCommand.get()
    
    def isQueueMandatory(self):
        return self.queueSystem.mandatory.get()
    
    def getSubmitTemplate(self):
        return self.queueSystem.getSubmitTemplate()
    
    def getQueuesDefault(self):
        return self.queueSystem.queuesDefault
    
    def getMpiCommand(self):
        return self.mpiCommand.get()
    
    def getQueueSystem(self):
        return self.queueSystem
    
    def setLabel(self, label):
        self.label.set(label)
    
    def setHostName(self, hostName):
        self.hostName.set(hostName)
    
    def setUserName(self, userName):
        self.userName.set(userName)
    
    def setPassword(self, password):
        self.password.set(password)
    
    def setHostPath(self, hostPath):
        self.hostPath.set(hostPath)
    
    def setMpiCommand(self, mpiCommand):
        self.mpiCommand.set(mpiCommand)  
        
    def setQueueSystem(self, queueSystem):
        self.queueSystem = queueSystem
        
    def getScipionHome(self):
        """ Return the path where Scipion is installed in 
        the host. This is useful when launching remote jobs.
        """ 
        return self.scipionHome.get()
    
    def setScipionHome(self, newScipionHome):
        self.scipionHome.set(newScipionHome)
        
    def getScipionConfig(self):
        """ From which file to read the configuration file in 
        this hosts. Useful for remote jobs.
        """
        return self.scipionConfig.get()
    
    def setScipionConfig(self, newConfig):
        self.scipionConfig.set(newConfig)
        
    def getAddress(self):
        return self.address.get()
    
    def setAddress(self, newAddress):
        return self.address.set(newAddress)


class QueueSystemConfig(OrderedObject):
    def __init__(self, **kwargs):
        OrderedObject.__init__(self, **kwargs)
        self.name = String()
        # Number of cores from which the queue is mandatory
        # 0 means no mandatory at all
        # 1 will force to launch all jobs through the queue
        self.mandatory = Integer()
        self.queues = None # List for queue configurations
        self.submitCommand = String()
        # Allow to change the prefix of submission scripts
        # we used by default the ID.job, but in some clusters
        # the job script should start by a letter
        self.submitPrefix = String()
        self.checkCommand = String()
        self.cancelCommand = String()
        self.submitTemplate = String()
        
    def hasName(self):
        return self.name.hasValue()
    
    def hasValue(self):
        return self.hasName() and len(self.queues)
    
    def getName(self):
        return self.name.get()
    
    def getMandatory(self):
        return self.mandatory.get()
    
    def getSubmitTemplate(self):
        return self.submitTemplate.get()
    
    def getSubmitCommand(self):
        return self.submitCommand.get()
    
    def getCheckCommand(self):
        return self.checkCommand.get()
    
    def getCancelCommand(self):
        return self.cancelCommand.get()
    
    def getQueues(self):
        return self.queues
    
    def setName(self, name):
        self.name.set(name)
    
    def setMandatory(self, mandatory):
        # This condition is to be backward compatible
        # when mandatory was a boolean
        # now it should use the number of CPU
        # that should force to use the queue
        if mandatory in ['False', 'false']:
            mandatory = 0
        elif mandatory in ['True', 'true']:
            mandatory = 1
            
        self.mandatory.set(mandatory)
    
    def setSubmitTemplate(self, submitTemplate):
        self.submitTemplate.set(submitTemplate)
    
    def setSubmitCommand(self, submitCommand):
        self.submitCommand.set(submitCommand)
    
    def setCheckCommand(self, checkCommand):
        self.checkCommand.set(checkCommand)
    
    def setCancelCommand(self, cancelCommand):
        self.cancelCommand.set(cancelCommand)
    
    def setQueues(self, queues):
        self.queues = queues
        
    def getQueueConfig(self, objId):
        if objId is not None and self.queues is not None:
            for queueConfig in self.queues:
                if objId == queueConfig.getObjId():
                    return queueConfig
        return None
        
        
#TODO: maybe deprecated
class QueueConfig(OrderedObject):
    def __init__(self, **kwargs):
        OrderedObject.__init__(self, **kwargs)
        self.name = String('default')
        self.maxCores = Integer()
        self.allowMPI = Boolean()
        self.allowThreads = Boolean()
        self.maxHours = Integer()
        
    def getName(self):
        return self.name.get()
    
    def getMaxCores(self):
        return self.maxCores.get()
    
    def getAllowMPI(self):
        return self.allowMPI.get()
    
    def getAllowThreads(self):
        return self.allowThreads.get()
    
    def getMaxHours(self):
        return self.maxHours.get()
    
    def setName(self, name):
        self.name.set(name)
    
    def setMaxCores(self, maxCores):
        self.maxCores.set(maxCores)
        
    def setAllowMPI(self, allowMPI):
        self.allowMPI.set(allowMPI)
    
    def setAllowThreads(self, allowThreads):
        self.allowThreads.set(allowThreads)
    
    def setMaxHours(self, maxHours):
        self.maxHours.set(maxHours)
   
