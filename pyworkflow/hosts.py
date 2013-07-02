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
from pyworkflow.object import *
from pyworkflow.mapper import SqliteMapper, XmlMapper


# class ExecutionHostConfig(OrderedObject):
#     """ Main store the configuration for execution hosts. """
#     
#     def __init__(self, **args):
#         OrderedObject.__init__(self, **args)
#         self.label = String()
#         self.hostName = String()
#         self.userName = String()
#         self.password = String()
#         self.hostPath = String()
#         self.mpiCommand = String()
#         self.queueSystem = QueueSystemConfig()
#     
#     def getLabel(self):
#         return self.label.get()
#     
#     def getHostName(self):
#         return self.hostName.get()
#     
#     def getUserName(self):
#         return self.userName.get()
#     
#     def getPassword(self):
#         return self.password.get()
#     
#     def getHostPath(self):
#         return self.hostPath.get()
#     
#     def getSubmitCommand(self):
#         return self.queueSystem.submitCommand.get()
#     
#     def isQueueMandatory(self):
#         return self.queueSystem.mandatory.get()
#     
#     def setLabel(self, label):
#         return self.label.set(label)
#     
#     def setHostName(self, hostName):
#         return self.hostName.set(hostName)
#     
#     def setUserName(self, userName):
#         return self.userName.set(userName)
#     
#     def setPassword(self, password):
#         return self.password.set(password)
#     
#     def setHostPath(self, hostPath):
#         return self.hostPath.set(hostPath)

        
class QueueSystemConfig(OrderedObject):
    def __init__(self, **args):
        OrderedObject.__init__(self, **args) 
        self.name = String()
        self.mandatory = Boolean()
        self.queues = List() # List for queue configurations
        self.submitCommand = String()
        self.checkCommand = String()
        self.cancelCommand = String()
        self.submitTemplate = String()
        
    def hasValue(self):
        return self.name.hasValue() and len(self.queues)
        
        
class QueueConfig(OrderedObject):
    def __init__(self, **args):
        OrderedObject.__init__(self, **args) 
        self.name = String('default')
        self.maxCores = Integer()
        self.allowMPI = Boolean()
        self.allowThreads = Boolean()
        self.maxHours = Integer()
        
    
# class ExecutionHostMapper(XmlMapper):
#     def __init__(self, filename, dictClasses=None, **args):
#         if dictClasses is None:
#             dictClasses = globals()
#         XmlMapper.__init__(self, filename, dictClasses, **args)
#         #self.setClassTag('ExecutionHostConfig.QueueSystemConfig', 'class_name')
#         self.setClassTag('ExecutionHostConfig.String', 'attribute')
#         self.setClassTag('ExecutionHostConfig.hostName', 'name_only')
#         self.setClassTag('ExecutionHostConfig.userName', 'name_only')
#         self.setClassTag('ExecutionHostConfig.password', 'name_only')
#         self.setClassTag('ExecutionHostConfig.hostPath', 'name_only')
#         self.setClassTag('ExecutionHostConfig.mpiCommand', 'name_only')
#         self.setClassTag('QueueSystemConfig.name', 'attribute')
#         self.setClassTag('QueueSystemConfig.mandatory', 'attribute')  
#         self.setClassTag('QueueConfig.ALL', 'attribute')   
#         self.setClassTag('List.QueueConfig', 'class_only')           
#         #self.setClassTag('ProtocolConfig.ProtocolConfig', 'class_only')
#     
#     def selectByLabel(self, objLabel):
#         hostsList = self.selectAll()
#         for host in hostsList:
#             if host.label == objLabel:
#                 return host
#         return None
    
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
    
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
        self.label = String()
        self.hostName = String()
        self.userName = String()
        self.password = String()
        self.hostPath = String()
        self.mpiCommand = String()
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
    
    def isQueueMandatory(self):
        return self.queueSystem.mandatory.get()
    
    def getMpiCommand(self):
        return self.mpiCommand.get()
    
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
