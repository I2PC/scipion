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
import re
from os.path import abspath, split

from pyworkflow.em import *
from pyworkflow.apps.config import *
from pyworkflow.protocol import *
from pyworkflow.mapper import SqliteMapper
from pyworkflow.utils import cleanPath, makePath, makeFilePath, join, exists, runJob, copyFile, timeit
from pyworkflow.utils.graph import Graph
from pyworkflow.hosts import HostMapper, HostConfig
import pyworkflow.protocol.launch as jobs
from pyworkflow.em.constants import RELATION_TRANSFORM, RELATION_SOURCE

PROJECT_DBNAME = 'project.sqlite'
PROJECT_LOGS = 'Logs'
PROJECT_RUNS = 'Runs'
PROJECT_TMP = 'Tmp'
PROJECT_SETTINGS = 'settings.sqlite'

# Regex to get numbering suffix and automatically propose runName
REGEX_NUMBER_ENDING = re.compile('(?P<prefix>.+\D)(?P<number>\d*)\s*$')


        
class Project(object):
    """This class will handle all information 
    related with a Project"""
    def __init__(self, path):
        """Create a project associated with a given path"""
        # To create a Project, a path is required
        self.name = path
        self.path = abspath(path)
        self.pathList = [] # Store all related paths
        self.dbPath = self.addPath(PROJECT_DBNAME)
        self.logsPath = self.addPath(PROJECT_LOGS)
        self.runsPath = self.addPath(PROJECT_RUNS)
        self.tmpPath = self.addPath(PROJECT_TMP)
        self.settingsPath = self.addPath(PROJECT_SETTINGS)
        self.runs = None
        self._runsGraph = None
        self._transformGraph = None
        self._sourceGraph = None
        
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
    
    def getDbPath(self):
        """ Return the path to the sqlite db. """
        return self.dbPath
    
    def getName(self):
        return self.name
    
    def getTmpPath(self, *paths):
        return self.getPath(PROJECT_TMP, *paths)
    
    def getSettings(self):
        return self.settings
    
    def saveSettings(self):
        # Read only mode
        if not isReadOnly():
            self.settings.write()
            
    def load(self):
        """Load project data and settings
        from the project dir."""
        
        if not exists(self.path):
            raise Exception("Cannot load project, path doesn't exist: %s" % self.path)
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
        self._cleanData()
        print abspath(self.dbPath)
        # Create db throught the mapper
        self.mapper = SqliteMapper(self.dbPath, globals())
        self.mapper.commit()
        # Write settings to disk
        self.settings = defaultSettings
        
        # Read only mode
        if not isReadOnly():
            
            self.settings.write(self.settingsPath)
            # Create other paths inside project
            for p in self.pathList:
                if '.' in p:
                    makeFilePath(p)
                else:
                    makePath(p)
        
    def _cleanData(self):
        """Clean all project data"""
        # Read only mode
        if not isReadOnly():
            cleanPath(*self.pathList)      
                
    def launchProtocol(self, protocol, wait=False):
        """ In this function the action of launching a protocol
        will be initiated. Actions done here are:
        1. Store the protocol and assign name and working dir
        2. Create the working dir and also the protocol independent db
        3. Call the launch method in protocol.job to handle submition: mpi, thread, queue,
        and also take care if the execution is remotely."""
        
        #if protocol.getStatus() != STATUS_INTERACTIVE:
        if not protocol.isInteractive.get():
            self._checkModificationAllowed([protocol], 'Cannot RE-LAUNCH protocol')
        
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
        
    def _updateProtocol(self, protocol, tries=0):
        # Read only mode
        if not isReadOnly():
            try:
                # FIXME: this will not work for a real remote host
                jobId = protocol.getJobId() # Preserve the jobId before copy
                dbPath = self.getPath(protocol.getDbPath())
                #join(protocol.getHostConfig().getHostPath(), protocol.getDbPath())
                prot2 = getProtocolFromDb(dbPath, protocol.getObjId())
                # Copy is only working for db restored objects
                protocol.setMapper(self.mapper)
                protocol.copy(prot2, copyId=False)
                # Restore jobId
                protocol.setJobId(jobId)
                self.mapper.store(protocol)
            except Exception, ex:
                print "Error trying to update protocol: %s(jobId=%s)\n ERROR: %s, tries=%d" % (protocol.getName(), jobId, ex, tries)
                if tries == 2: # 3 tries have been failed
                    import traceback
                    traceback.print_exc()
                    # If any problem happens, the protocol will be marked wih a status fail
                    protocol.setFailed(str(ex))
                    self.mapper.store(protocol)
                else:
                    time.sleep(1)
                    self._updateProtocol(protocol, tries+1)
        
    def stopProtocol(self, protocol):
        """ Stop a running protocol """
        try:
            jobs.stop(protocol)
        except Exception:
            raise
        finally:
            protocol.setAborted()
            self._storeProtocol(protocol)
        
    def continueProtocol(self, protocol):
        """ This function should be called 
        to mark a protocol that have an interactive step
        waiting for approval that can continue
        """
        protocol.continueFromInteractive()
        self.launchProtocol(protocol)
        
    def __protocolInList(self, prot, protocols):
        """ Check if a protocol is in a list comparing the ids. """
        for p in protocols:
            if p.getObjId() == prot.getObjId():
                return True
        return False
    
    def __validDependency(self, prot, child, protocols):
        """ Check if the given child is a true dependency of the protocol
        in order to avoid any modification.
        """
        return (not self.__protocolInList(child, protocols) and
                not child.isSaved()) 
        
    
    def _checkProtocolsDependencies(self, protocols, msg):
        """ Check if the protocols have depencies.
        This method is used before delete or save protocols to be sure
        it is not referenced from other runs. (an Exception is raised)
        Params:
             protocols: protocol list to be analyzed.
             msg: String message to be prefixed to Exception error.
        """
        # Check if the protocol have any dependencies
        error = ''
        for prot in protocols:
            node = self.getRunsGraph().getNode(prot.strId())
            if node:
                childs = [node.run for node in node.getChilds() if self.__validDependency(prot, node.run, protocols)]
                if childs:
                    deps = [' ' + c.getRunName() for c in childs]
                    error += '\n *%s* is referenced from:\n   - ' % prot.getRunName()
                    error += '\n   - '.join(deps) 
        if error:
            raise Exception(msg + error)
        
    def _checkModificationAllowed(self, protocols, msg):
        """ Check if any modification operation is allowed for
        this group of protocols. 
        """
        if isReadOnly():
            raise Exception(msg + " Running in READ-ONLY mode.")
        
        self._checkProtocolsDependencies(protocols, msg)        
        
    def deleteProtocol(self, *protocols):
        self._checkModificationAllowed(protocols, 'Cannot DELETE protocols')
        
        for prot in protocols:
            self.mapper.delete(prot) # Delete from database
            wd = prot.workingDir.get()
            
            if wd.startswith(PROJECT_RUNS):
                cleanPath(wd)
            else:
                print "Error path: ", wd 
      
        self.mapper.commit()     
        
    def __setProtocolLabel(self, newProt):
        """ Set a readable label to a newly created protocol.
        We will try to find another existing protocol of the 
        same class and set the label from it.
        """
        prevLabel = None # protocol with same class as newProt
        newProtClass = newProt.getClass()
        
        for prot in self.getRuns(iterate=True, refresh=False):
            if newProtClass == prot.getClass():
                prevLabel = prot.getObjLabel().strip()
                
        if prevLabel: 
            numberSuffix = 2
            m = REGEX_NUMBER_ENDING.match(prevLabel)
            if m and m.groupdict()['number']:
                numberSuffix = int(m.groupdict()['number']) + 1
                prevLabel = m.groupdict()['prefix']
            protLabel =  prevLabel + ' %s' % numberSuffix
        else:
            protLabel = newProt.getClassLabel()
            
        newProt.setObjLabel(protLabel)
        
    def newProtocol(self, protocolClass, **kwargs):
        """ Create a new protocol from a given class. """
        newProt = protocolClass(project=self, **kwargs)
        self.__setProtocolLabel(newProt)
        
        newProt.setMapper(self.mapper)
        newProt.setProject(self)
        
        return newProt
                      
    def __getIOMatches(self, node, childNode):
        """ Check if some output of node is used as input in childNode.
        Return the list of attribute names that matches.
        Used from self.copyProtocol
        """
        matches = []
        for oKey, oAttr in node.run.iterOutputAttributes(EMObject):
            for iKey, iAttr in childNode.run.iterInputAttributes():
                if oAttr is iAttr.get():
                    matches.append((oKey, iKey))
        
        return matches                    
        
    def copyProtocol(self, protocol):
        """ Make a copy of the protocol, return a new one with copied values. """
        result = None
        
        if isinstance(protocol, Protocol):
            newProt = self.newProtocol(protocol.getClass())
            newProt.copyDefinitionAttributes(protocol)
            result = newProt
    
        elif isinstance(protocol, list):
            # Handle the copy of a list of protocols
            # for this case we need to update the references of input/outputs
            newDict = {}
                        
            for prot in protocol:
                newProt = self.newProtocol(prot.getClass())
                newProt.copyDefinitionAttributes(prot)
                newDict[prot.getObjId()] = newProt
                self.saveProtocol(newProt)
                         
            g = self.getRunsGraph(refresh=False)
            
            for prot in protocol:
                node = g.getNode(prot.strId())
                newProt = newDict[prot.getObjId()]
                
                for childNode in node.getChilds():
                    newChildProt = newDict.get(childNode.run.getObjId(), None)
                    
                    if newChildProt:
                        # Get the matches between outputs/inputs of node and childNode
                        matches = self.__getIOMatches(node, childNode)
                        #print "%s -> %s, matches: %s" % (prot.getRunName(), childNode.runmatches
                        # For each match, set the pointer and the extend attribute
                        # to reproduce the dependencies in the new workflow
                        for oKey, iKey in matches:
                            childPointer = getattr(newChildProt, iKey)
                            childPointer.set(newProt)
                            childPointer.setExtendedAttribute(oKey)
                        self.mapper.store(newChildProt)                   

            self.mapper.commit()
        else:
            raise Exception("Project.copyProtocol: invalid input protocol type '%s'." % type(protocol))
    
        return result
    
    def saveProtocol(self, protocol):
        self._checkModificationAllowed([protocol], 'Cannot SAVE protocol')
        
        protocol.setStatus(STATUS_SAVED)
        if protocol.hasObjId():
            self._storeProtocol(protocol)
        else:
            self._setupProtocol(protocol)
                
    def getProtocol(self, protId):
        return self.mapper.selectById(protId)
        
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
        # Read only mode
        if not isReadOnly():
            self.mapper.store(protocol)
            self.mapper.commit()
    
    def _setupProtocol(self, protocol):
        """Insert a new protocol instance in the database"""
        
        # Read only mode
        if not isReadOnly():
        
            self._storeProtocol(protocol) # Store first to get a proper id
            # Set important properties of the protocol
            name = protocol.getClassName() + protocol.strId()
            protocol.setProject(self)
            print protocol.strId(), protocol.getProject().getName()
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
                r.setProject(self)
                r.setMapper(self.mapper)
                # Update nodes that are running and are not invoked by other protocols
                if r.isActive():
                    if not r.isChild():
                        self._updateProtocol(r)
            self.mapper.commit()
        
        return self.runs
    
    def iterSubclasses(self, classesName, objectFilter=None):
        """ Retrieve all objects from the project that are instances
            of any of the classes in classesName list.
        Params: 
            classesName: String with commas separated values of classes name. 
            objectFilter: a filter function to discard some of the retrieved objects."""
        for objClass in classesName.split(","):
            for obj in self.mapper.selectByClass(objClass.strip(), iterate=True, objectFilter=objectFilter):
                yield obj
    
    def getRunsGraph(self, refresh=True):
        """ Build a graph taking into account the dependencies between
        different runs, ie. which outputs serves as inputs of other protocols. 
        """
        if refresh or self._runsGraph is None:
            outputDict = {} # Store the output dict
            runs = [r for r in self.getRuns(refresh=refresh) if not r.isChild()]
            g = Graph(rootName='PROJECT')
            
            for r in runs:
                n = g.createNode(r.strId())
                n.run = r
                n.label = r.getRunName()
                outputDict[r.getName()] = n
                for _, attr in r.iterOutputAttributes(EMObject):
                    outputDict[attr.getName()] = n # mark this output as produced by r
                
            def _checkInputAttr(node, pointed):
                """ Check if an attr is registered as output"""
                if pointed:
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
                        pointed = attr.getObjValue()
                        # Only checking pointed object and its parent, if more levels
                        # we need to go up to get the correct dependencies
                        (_checkInputAttr(node, pointed) or 
                         _checkInputAttr(node, self.mapper.getParent(pointed))
                        )
            rootNode = g.getRoot()
            rootNode.run = None
            rootNode.label = "PROJECT"
            
            for n in g.getNodes():
                if n.isRoot() and not n is rootNode:
                    rootNode.addChild(n)
            self._runsGraph = g
            
        return self._runsGraph
    
    def _getRelationGraph(self, relation=RELATION_SOURCE, refresh=False):
        """ Retrieve objects produced as outputs and
        make a graph taking into account the SOURCE relation. """
        relations = self.mapper.getRelationsByName(relation)
        g = Graph(rootName='PROJECT')
        root = g.getRoot()
        root.object = None
        runs = self.getRuns(refresh=refresh)
        
        for r in runs:
            for _, attr in r.iterOutputAttributes(EMObject):
                node = g.createNode(attr.strId(), attr.getNameId())
                node.object = attr                
        
        for rel in relations:
            pid = str(rel['object_parent_id'])
            parent = g.getNode(pid)
            if not parent:
                print "error, parent none: ", pid
            else:
                cid = str(rel['object_child_id'])
                child = g.getNode(cid)
                if not child:
                    print "error, child none: ", cid, " label: ", 
                    print "   parent: ", pid 
                else:
                    parent.addChild(child)
            
        for n in g.getNodes():
            if n.isRoot() and not n is root:
                root.addChild(n)
            
        return g
               
    def getTransformGraph(self, refresh=False):
        """ Get the graph from the TRASNFORM relation. """
        if refresh or not self._transformGraph:
            self._transformGraph = self._getRelationGraph(RELATION_TRANSFORM, refresh)
            
        return self._transformGraph
            
    def getSourceGraph(self, refresh=False):
        """ Get the graph from the SOURCE relation. """
        if refresh or not self._sourceGraph:
            self._sourceGraph = self._getRelationGraph(RELATION_SOURCE, refresh)
            
        return self._sourceGraph
                
    def getRelatedObjects(self, relation, obj, direction=RELATION_CHILDS):
        """ Get all objects related to obj by a give relation.
        Params:
            relation: the relation name to search for.
            obj: object from which the relation will be search,
                actually not only this, but all other objects connected
                to this one by the RELATION_TRANSFORM.
            direction: this say if search for childs or parents in the relation.
        """
        graph = self.getTransformGraph()
        relations = self.mapper.getRelationsByName(relation)
        connection = self._getConnectedObjects(obj, graph)
        objects = []
        
        for rel in relations:
            pid = str(rel['object_parent_id'])
            parent = connection.get(pid, None)
            if parent:
                cid = str(rel['object_child_id'])
                child = graph.getNode(cid).object
                objects.append(child)
                
        return objects
    
    def _getConnectedObjects(self, obj, graph):
        """ Give a TRANSFORM graph, return the elements that
        are connected, either childs, ancestors or siblings. 
        """
        n = graph.getNode(obj.strId())
        # Get the oldest ancestor of a node, before 
        # reaching the root node
        while not n.getParent().isRoot():
            n = n.getParent()
            
        connection = {}
        for node in n.iterChilds():
            connection[node.getName()] = node.object
        
        return connection
            

def isReadOnly():
    """ Auxiliar method to keep a read-only mode for the environment. """
    return 'SCIPION_READONLY' in os.environ
        
