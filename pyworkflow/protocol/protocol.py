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
This modules contains classes required for the workflow
execution and tracking like: Step and Protocol
"""

import os
import sys
import datetime as dt
import pickle
import json
from collections import OrderedDict
import traceback
import time

import pyworkflow as pw
from pyworkflow.object import *
from pyworkflow.utils import redStr, greenStr, magentaStr, envVarOn, runJob
from pyworkflow.utils.path import (makePath, join, missingPaths, cleanPath, cleanPattern,
                                   getFiles, exists, renderTextFile, copyFile)
from pyworkflow.utils.log import ScipionLogger
from executor import StepExecutor, ThreadStepExecutor, MPIStepExecutor
from constants import *
from params import Form
import scipion



class Step(OrderedObject):
    """ Basic execution unit.
    It should defines its Input, Output
    and define a run method.
    """
    def __init__(self, **kwargs):
        OrderedObject.__init__(self, **kwargs)
        self._prerequisites = CsvList() # which steps needs to be done first
        self.status = String()
        self.initTime = String()
        self.endTime = String()
        self._error = String()
        self.interactive = Boolean(False)
        self._resultFiles = String()
        self._index = None

    def getIndex(self):
        return self._index
    
    def setIndex(self, newIndex):
        self._index = newIndex
        
    def getPrerequisites(self):
        return self._prerequisites
    
    def addPrerequisites(self, *newPrerequisites):
        for p in newPrerequisites:
            self._prerequisites.append(p)

    def _preconditions(self):
        """ Check if the necessary conditions to
        step execution are met""" 
        return self._validate() == []
    
    def _postconditions(self):
        """ Check if the step have done well its task
        and accomplish its results""" 
        return True
    
    def _run(self):
        """ This is the function that will do the real job.
        It should be override by sub-classes.""" 
        pass
    
    def setRunning(self):
        """ The the state as STATE_RUNNING and 
        set the init and end times.
        """
        self.initTime.set(dt.datetime.now())
        self.endTime.set(None)
        self.status.set(STATUS_RUNNING)
        self._error.set(None) # Clean previous error message
        
    def getError(self):
        return self._error
    
    def getErrorMessage(self):
        return self.getError().get('')
        
    def setFailed(self, msg):
        """ Set the run failed and store an error message. """
        self.endTime.set(dt.datetime.now())
        self._error.set(msg)
        self.status.set(STATUS_FAILED)
        
    def setAborted(self):
        """ Set the status to aborted and updated the endTime. """
        self.endTime.set(dt.datetime.now())
        self._error.set("Aborted by user.")
        self.status.set(STATUS_ABORTED)

    def getStatus(self):
        return self.status.get(STATUS_NEW)
    
    def getElapsedTime(self, default=dt.timedelta()):
        """ Return the time that took to run 
        (or the actual running time if still is running )
        """
        elapsed = default

        if self.initTime.hasValue():
            t1 = self.initTime.datetime()

            if self.endTime.hasValue():
                t2 = self.endTime.datetime()
            else:
                t2 = dt.datetime.now()

            elapsed = t2 - t1
        
        return elapsed
    
    def setStatus(self, value):
        return self.status.set(value)
    
    def setInteractive(self, value):
        return self.interactive.set(value)
    
    def isActive(self):
        return self.getStatus() in ACTIVE_STATUS
    
    def isFinished(self):
        return self.getStatus() == STATUS_FINISHED
    
    def isRunning(self):
        return self.getStatus() == STATUS_RUNNING
    
    def isFailed(self):
        return self.getStatus() == STATUS_FAILED
    
    def isSaved(self):
        return self.getStatus() == STATUS_SAVED
    
    def isAborted(self):
        return self.getStatus() == STATUS_ABORTED

    def isLaunched(self):
        return self.getStatus() == STATUS_LAUNCHED

    def isInteractive(self):
        return self.interactive.get()
    
    def isWaiting(self):
        return self.getStatus() == STATUS_WAITING

    def run(self):
        """ Do the job of this step"""
        self.setRunning() 
        try:
            self._run()
            self.endTime.set(dt.datetime.now())
            if self.status.get() == STATUS_RUNNING:
                if self.isInteractive():
                    # If the Step is interactive, after run
                    # it will be waiting for use to mark it as DONE
                    status = STATUS_INTERACTIVE
                else:
                    status = STATUS_FINISHED
                self.status.set(status)
        except Exception, e:
            self.setFailed(str(e))
            import traceback
            traceback.print_exc()            
            #raise #only in development
        # finally:
        #     self.endTime.set(dt.datetime.now())


class FunctionStep(Step):
    """ This is a Step wrapper around a normal function
    This class will ease the insertion of Protocol function steps
    through the function _insertFunctionStep"""
    def __init__(self, func=None, funcName=None, *funcArgs, **kwargs):
        """ 
         Params:
            func: the function that will be executed.
            funcName: the name assigned to that function (will be stored)
            *funcArgs: argument list passed to the function (serialized and stored)
            **kwargs: extra parameters.
        """ 
        Step.__init__(self)
        self._func = func # Function should be set before run
        self._args = funcArgs
        self.funcName = String(funcName)
        self.argsStr = String(pickle.dumps(funcArgs))
        self.setInteractive(kwargs.get('interactive', False))
        if kwargs.get('wait', False):
            self.setStatus(STATUS_WAITING)
        
    def _runFunc(self):
        """ Return the possible result files after running the function. """
        return self._func(*self._args)

    def _run(self):
        """ Run the function and check the result files if any. """
        resultFiles = self._runFunc() 
        if isinstance(resultFiles, basestring):
            resultFiles = [resultFiles]
        if resultFiles and len(resultFiles):
            missingFiles = missingPaths(*resultFiles)
            if len(missingFiles):
                raise Exception('Missing filePaths: ' + ' '.join(missingFiles))
            self._resultFiles.set(pickle.dumps(resultFiles))
    
    def _postconditions(self):
        """ This type of Step, will simply check
        as postconditions that the result filePaths exists""" 
        if not self._resultFiles.hasValue():
            return True
        filePaths = pickle.loads(self._resultFiles.get())

        return len(missingPaths(*filePaths)) == 0
    
    def __eq__(self, other):
        """ Compare with other FunctionStep""" 
        return (self.funcName == other.funcName and
                self.argsStr == other.argsStr)
        
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __str__(self):
        return self.funcName.get()
        
            
class RunJobStep(FunctionStep):
    """ This Step will wrapper the commonly used function runJob
    for launching specific programs with some parameters.
    The runJob function should be provided by the protocol
    when inserting a new RunJobStep"""
     
    def __init__(self, runJobFunc=None, programName=None, arguments=None, resultFiles=[], **kwargs):
        FunctionStep.__init__(self, runJobFunc, 'runJob', programName, arguments)
        # Number of mpi and threads used to run the program
        self.__runJob = runJobFunc # Store the current function to run the job
        self.mpi = 1
        self.threads = 1
        
    def _runFunc(self):
        """ Wrap around runJob function""" 
        # We know that:
        #  _func: is the runJob function
        #  _args[0]: is the program name
        #  _args[1]: is the argumets to the program
        return self._func(None, self._args[0], self._args[1], 
               numberOfMpi=self.mpi, numberOfThreads=self.threads)
        #TODO: Add the option to return resultFiles
    
    def __str__(self):
        return self._args[0] # return program name         
                
                
class StepSet(Set):
    """ Special type of Set for storing steps. """
    def __init__(self, filename=None, prefix='', 
                 mapperClass=None, **kwargs):
        Set.__init__(self, filename, prefix, mapperClass, classesDict=globals(), **kwargs)
        
                
class Protocol(Step):
    """ The Protocol is a higher type of Step.
    It also have the inputs, outputs and other Steps properties,
    but contains a list of steps that are executed
    """
    def __init__(self, **kwargs):
        Step.__init__(self, **kwargs)        
        self._steps = [] # List of steps that will be executed
        self.workingDir = String(kwargs.get('workingDir', '.')) # All generated filePaths should be inside workingDir
        self.mapper = kwargs.get('mapper', None)
        self._inputs = []
        self._outputs = CsvList()
        # Expert level
        self.expertLevel = Integer(kwargs.get('expertLevel', LEVEL_NORMAL))#needs to be defined before parsing params 
        self._definition = Form(self)
        self._defineParams(self._definition)
        self._createVarsFromDefinition(**kwargs)
        self.__stdOut = None
        self.__stdErr = None
        self.__fOut = None
        self.__fErr = None
        self._log = None
        self._buffer = ''  # text buffer for reading log files
        # Project to which the protocol belongs
        self.__project = kwargs.get('project', None)
        # Filename templates dict that will be used by _getFileName
        self.__filenamesDict = {}
        
        # For non-parallel protocols mpi=1 and threads=1
        self.allowMpi = hasattr(self, 'numberOfMpi')
        if not self.allowMpi:
            self.numberOfMpi = Integer(1)
        
        self.allowThreads = hasattr(self, 'numberOfThreads')
        
        if not self.allowThreads:
            self.numberOfThreads = Integer(1)
        
        # Check if MPI or threads are passed in **kwargs, mainly used in tests
        if 'numberOfMpi' in kwargs:
            self.numberOfMpi.set(kwargs.get('numberOfMpi'))
        
        if 'numberOfThreads' in kwargs:
            self.numberOfThreads.set(kwargs.get('numberOfThreads'))            
        
        if not hasattr(self, 'hostName'):
            self.hostName = String(kwargs.get('hostName', 'localhost'))
        
        # Maybe this property can be inferred from the 
        # prerequisites of steps, but is easier to keep it
        self.stepsExecutionMode = STEPS_SERIAL
        
        # Run mode
        self.runMode = Integer(kwargs.get('runMode', MODE_RESUME))
        # Use queue system?
        self._useQueue = Boolean(False)
        # Store a json string with queue name
        # and queue parameters (only meanful if _useQueue=True)
        self._queueParams = String()
        
        self._jobId = String() # Store queue job id
        self._pid = Integer()
        self._stepsExecutor = None
        self._stepsDone = Integer(0)
        self._numberOfSteps = Integer(0)
        # For visualization
        self.allowHeader = Boolean(True)    
        # Create an String variable to allow some protocol to precompute
        # the summary message
        self.summaryVar = String()
        self.methodsVar = String()    
        
    def _storeAttributes(self, attrList, attrDict):
        """ Store all attributes in attrDict as 
        attributes of self, also store the key in attrList.
        """
        for key, value in attrDict.iteritems():
            if key not in attrList:
                attrList.append(key)
            setattr(self, key, value)
        
    def _defineInputs(self, **kwargs):
        """ This function should be used to define
        those attributes considered as Input.
        """ 
        self._storeAttributes(self._inputs, kwargs)
        
    def _defineOutputs(self, **kwargs):
        """ This function should be used to specify
        expected outputs.
        """ 
        for k, v in kwargs.iteritems():
            if hasattr(self, k):
                self._deleteChild(k, v)
            self._insertChild(k, v)
            
        self._storeAttributes(self._outputs, kwargs)
        
    def _updateOutputSet(self, outputName, outputSet, 
                         state=Set.STREAM_OPEN):
        """ Use this function when updating an Stream output set.
        """
        self.__tryUpdateOuputSet(outputName, outputSet, state)

    def __tryUpdateOuputSet(self, outputName, outputSet,
                         state=Set.STREAM_OPEN, tries=1):
        try:
            # Update the set with the streamState value (either OPEN or CLOSED)
            outputSet.setStreamState(state)

            if self.hasAttribute(outputName):
                outputSet.write() # Write to commit changes
                outputAttr = getattr(self, outputName)
                # Copy the properties to the object contained in the protcol
                outputAttr.copy(outputSet, copyId=False)
                # Persist changes
                self._store(outputAttr)
            else:
                # Here the defineOutputs function will call the write() method
                self._defineOutputs(**{outputName: outputSet})
                self._store(outputSet)
            # Close set databaset to avoid locking it
            outputSet.close()

        except Exception as ex:
            print("Error trying to update output of protocol, tries=%d" % tries)

            if tries > 3:
                raise ex
            else:
                time.sleep(tries)
                self.__tryUpdateOuputSet(outputName, outputSet, state, tries+1)
        
    def getProject(self):
        return self.__project
    
    def setProject(self, project):
        self.__project = project
        
    @staticmethod
    def hasDefinition(cls):
        """ Check if the protocol has some definition.
        This can help to detect "abstract" protocol that 
        only serve as base for other, not to be instantiated. 
        """
        return hasattr(cls, '_definition')
    
    def getDefinition(self):
        """ Access the protocol definition. """
        return self._definition
    
    def getParam(self, paramName):
        """ Return a _definition param give its name. """
        return self._definition.getParam(paramName)
    
    def getEnumText(self, paramName):
        """ This function will retrieve the text value
        of an enum parameter in the definition, taking the actual value in the protocol.
        Params:
            paramName: the name of the enum param.
        Returns:
            the string value corresponding to the enum choice.
        """
        index = getattr(self, paramName).get() # self.getAttributeValue(paramName)
        return self.getParam(paramName).choices[index]
    
    def evalParamCondition(self, paramName):
        """ Eval if the condition of paramName in _definition
        is satified with the current values of the protocol attributes. 
        """
        return self._definition.evalParamCondition(paramName)
    
    def evalExpertLevel(self, paramName):
        """ Return the expert level evaluation for a param with the given name. """
        return self.evalParamExpertLevel(self.getParam(paramName))
    
    def evalParamExpertLevel(self, param):
        """ Return True if the param has an expert level is less than 
        the one for the whole protocol. 
        """
        return param.expertLevel.get() <= self.expertLevel.get()
    
    def iterDefinitionAttributes(self):
        """ Iterate over all the attributes from definition. """
        for paramName, _ in self._definition.iterParams():
            yield paramName, getattr(self, paramName)
            
    def getDefinitionDict(self):
        """ Similar to getObjDict, but only for those 
        params that are in the form.
        This function is used for export protocols as json text file.
        """
        d = OrderedDict()
        d['object.className'] = self.getClassName()
        d['object.id'] = self.strId()
        d['object.label'] = self.getObjLabel()
        d['object.comment'] = self.getObjComment()
         
        od = self.getObjDict()
        
        for attrName in od:
            if self.getParam(attrName) is not None:
                d[attrName] = od[attrName]
                
        return d
        
    def iterDefinitionSections(self):
        """ Iterate over all the section of the definition. """
        for section in self._definition.iterSections():
            yield section
            
    def iterInputAttributes(self):
        """ Iterate over the main input parameters
        of this protocol. Now the input are assumed to be these attribute
        which are pointers and have no condition.
        """
        for key, attr in self.getAttributes():
            if not isinstance(attr, Object):
                raise Exception('Attribute %s have been overwritten to type %s ' % (key, type(attr)))
            if isinstance(attr, PointerList) and attr.hasValue():
                for item in attr:
                    # the same key is returned for all items inside the
                    # PointerList, this is used in viewprotocols.py
                    # to group them inside the same tree element
                    yield key, item 
            elif attr.isPointer() and attr.hasValue():
                yield key, attr
                
    def iterOutputAttributes(self, outputClass):
        """ Iterate over the outputs produced by this protocol. """
        for key, attr in self.getAttributes():
            if isinstance(attr, outputClass):
                yield key, attr
                
    def iterOutputEM(self):
        """ Iterate over the outputs that are EM objects. """
        from pyworkflow.em.data import EMObject
        for paramName, attr in self.iterOutputAttributes(EMObject):
            yield paramName, attr

    def getOutputsSize(self):
        return sum(1 for _ in self.iterOutputEM())

    def getOutputFiles(self):
        """ Return the output files produced by this protocol.
        This can be used in web to download or in remote 
        executions to copy results back.
        """
        # By default return the output file of each output attribute
        s = set()
        
        for _, attr in self.iterOutputEM():
            s.update(attr.getFiles())
            
        return s
            
    def copyDefinitionAttributes(self, other):
        """ Copy definition attributes to other protocol. """
        for paramName, _ in self.iterDefinitionAttributes():
            self.copyAttributes(other, paramName)
        
    def _createVarsFromDefinition(self, **kwargs):
        """ This function will setup the protocol instance variables
        from the Protocol Class definition, taking into account
        the variable type and default values.
        """
        if hasattr(self, '_definition'):
            for paramName, param in self._definition.iterParams():
                # Create the var with value coming from kwargs or from 
                # the default param definition
                var = param.paramClass(value=kwargs.get(paramName, param.default.get()))
                setattr(self, paramName, var)
        else:
            print "FIXME: Protocol '%s' has not DEFINITION" % self.getClassName()
        
    def _getFileName(self, key, **kwargs):
        """ This function will retrieve filenames given a key and some
        keywords arguments. The __filenamesDict attribute should be 
        updated with templates that accept the given keys.
        """
        return self.__filenamesDict[key] % kwargs
    
    def _updateFilenamesDict(self, fnDict):
        """ Update the dictionary with templates that will be used
        by the _getFileName function.
        """
        self.__filenamesDict.update(fnDict)
    
    def _store(self, *objs):
        """ Stores objects of the protocol using the mapper.
        If not objects are passed, the whole protocol is stored.
        """
        if not self.mapper is None:
            if len(objs) == 0:
                self.mapper.store(self)
            else:
                for obj in objs:
                    self.mapper.store(obj)
            self.mapper.commit()
            
    def _insertChild(self, key, child):
        """ Insert a new child not stored previously.
        If stored previously, _store should be used.
        The child will be set as self.key attribute
        """
        try:
            setattr(self, key, child)
            if self.hasObjId():
                self.mapper.insertChild(self, key, child)
        except Exception, ex:
            print "Error with child '%s', value=%s, type=%s" % (key, child, type(child))
            raise ex
        
    def _deleteChild(self, key, child):
        """ Delete a child from the mapper. """
        self.mapper.delete(child)
        
    def _insertAllSteps(self):
        """ Define all the steps that will be executed. """
        pass
    
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        pass
    
    def __insertStep(self, step, **kwargs):
        """ Insert a new step in the list. 
        Params:
         **kwargs:
            prerequisites: a list with the steps index that need to be done 
                           previous than the current one."""
        prerequisites = kwargs.get('prerequisites', None)
        
        if prerequisites is None:
            if len(self._steps):
                step.addPrerequisites(len(self._steps)) # By default add the previous step as prerequisite
        else:
            step.addPrerequisites(*prerequisites)
                
        self._steps.append(step)
        # Setup and return step index
        step.setIndex(len(self._steps))
                
        return step.getIndex()
        
    def _getPath(self, *paths):
        """ Return a path inside the workingDir. """
        return join(self.workingDir.get(), *paths)

    def _getExtraPath(self, *paths):
        """ Return a path inside the extra folder. """
        return self._getPath("extra", *paths)
    
    def _getTmpPath(self, *paths):
        """ Return a path inside the tmp folder. """
        return self._getPath("tmp", *paths) 
    
    def _getLogsPath(self, *paths):
        return self._getPath("logs", *paths)   
    
    def _getRelPath(self, *path):
        """ Return a relative path from the workingDir. """
        return os.path.relpath(self._getPath(*path), self.workingDir.get())
    
    def _getRelPathExecutionDir(self, *path):
        """ Return a relative path from the projdir. """
        #TODO  must be a bettis
        return os.path.relpath(
                  self._getPath(*path),
                  os.path.dirname(os.path.dirname(self.workingDir.get()))
                              )   
    
    def _getBasePath(self, path):
        """ Take the basename of the path and get the path
        relative to working dir of the protocol. 
        """
        return self._getPath(os.path.basename(path))
        
    def _insertFunctionStep(self, funcName, *funcArgs, **kwargs):
        """ 
         Params:
           funcName: the string name of the function to be run in the Step.
           *funcArgs: the variable list of arguments to pass to the function.
           **kwargs: see __insertStep
        """
        # Get the function give its name
        func = getattr(self, funcName, None)
        
        # Ensure the protocol instance have it and is callable
        if not func:
            raise Exception("Protocol._insertFunctionStep: '%s' function is not member of the protocol" % funcName)
        if not callable(func):
            raise Exception("Protocol._insertFunctionStep: '%s' is not callable" % funcName)
        step = FunctionStep(func, funcName, *funcArgs, **kwargs)
        
        return self.__insertStep(step, **kwargs)
        
    def _insertRunJobStep(self, progName, progArguments, resultFiles=[], **kwargs):
        """ Insert an Step that will simple call runJob function
        **args: see __insertStep
        """
        return self._insertFunctionStep('runJob', progName, progArguments, **kwargs)
    
    def _insertCopyFileStep(self, sourceFile, targetFile, **kwargs):
        """ Shortcut function to insert an step for copying a file to a destiny. """
        step = FunctionStep(copyFile, 'copyFile', sourceFile, targetFile, **kwargs)
        return self.__insertStep(step, **kwargs)
            
    def _enterDir(self, path):
        """ Enter into a new directory path and store the current path.
        The current path will be used in _leaveDir, but nested _enterDir
        are not allowed since self._currentDir is overriden.
        """
        self._currentDir = os.getcwd()
        os.chdir(path)
        if self._log:
            self._log.info("Entered into dir: cd '%s'" % path)
        
    def _leaveDir(self): 
        """ This method should be called after a call to _enterDir
        to return to the previous location. 
        """
        os.chdir(self._currentDir)  
        if self._log:
            self._log.info("Returned to dir: cd '%s'" % self._currentDir)
                      
    def _enterWorkingDir(self):
        """ Change to the protocol working dir. """
        self._enterDir(self.workingDir.get())
        
    def _leaveWorkingDir(self):
        """ This funcion make sense to use in conjunction 
        with _enterWorkingDir to go back to execution path.
        """
        self._leaveDir()
        
    def continueFromInteractive(self):
        """ TODO: REMOVE this function.
        Check if there is an interactive step and set
        as finished, this is used now mainly in picking,
        but we should remove this since is weird for users.
        """
        if exists(self.getStepsFile()):
            stepsSet = StepSet(filename=self.getStepsFile())
            for step in stepsSet:
                if step.getStatus() == STATUS_INTERACTIVE:
                    step.setStatus(STATUS_FINISHED)
                    stepsSet.update(step)
                    break
            stepsSet.write()
            stepsSet.close() # Close the connection
    
    def loadSteps(self):
        """ Load the Steps stored in the steps.sqlite file.
        """
        prevSteps = []
        
        if exists(self.getStepsFile()):
            stepsSet = StepSet(filename=self.getStepsFile())
            for step in stepsSet:
                prevSteps.append(step.clone())
            stepsSet.close() # Close the connection
        return prevSteps
    
    def _insertPreviousSteps(self):
        """ Insert steps of previous execution. 
        It can be used to track previous steps done for
        protocol that allow some kind of continue (such as ctf estimation).
        """
        for step in self.loadSteps():
            self.__insertStep(step)
        
    def __findStartingStep(self):
        """ From a previous run, compare self._steps and self._prevSteps
        to find which steps we need to start at, skipping successful done 
        and not changed steps. Steps that needs to be done, will be deleted
        from the previous run storage.
        """
        if self.runMode == MODE_RESTART:
            self._prevSteps = []
            return 0
        
        self._prevSteps = self.loadSteps()
        
        n = min(len(self._steps), len(self._prevSteps))
        self.info("len(steps) " + str(len(self._steps)) + " len(prevSteps) " + str(len(self._prevSteps)))
        
        for i in range(n):
            newStep = self._steps[i]
            oldStep = self._prevSteps[i]
#             self.info(">>> i: %s" % i)
#             self.info(" oldStep.status: %s" % str(oldStep.status))
#             self.info(" oldStep!=newStep %s" % (oldStep != newStep))
#             self.info(" not post: %s" % (not oldStep._postconditions()))
            if (not oldStep.isFinished() or
                newStep != oldStep or 
                not oldStep._postconditions()):
#                 self.info(" returning i: %s" % i)
                return i
            newStep.copy(oldStep)
            
        return n
    
    def _storeSteps(self):
        """ Store the new steps list that can be retrieved 
        in further execution of this protocol.
        """
        stepsFn = self.getStepsFile()
           
        self._stepsSet = StepSet(filename=stepsFn)
        self._stepsSet.setStore(False)
        self._stepsSet.clear()
        
        for step in self._steps:
            step.cleanObjId()
            self.setInteractive(self.isInteractive() or step.isInteractive())
            self._stepsSet.append(step)

        self._stepsSet.write()
        
    def __updateStep(self, step):
        """ Store a given step and write changes. """
        self._stepsSet.update(step)
        self._stepsSet.write()
        
    def _stepStarted(self, step):
        """This function will be called whenever an step
        has started running.
        """
        self.info(magentaStr("STARTED") + ": %s, step %d" %
                  (step.funcName.get(), step._index))
        self.info("  %s" % step.initTime.datetime())
        self.__updateStep(step)
        
    def _stepFinished(self, step):
        """This function will be called whenever an step
        has finished its run.
        """
        doContinue = True
        if step.isInteractive():
            doContinue = False
        elif step.isFailed():
            doContinue = False
            errorMsg = redStr("Protocol failed: " + step.getErrorMessage())
            self.setFailed(errorMsg)
            self.error(errorMsg)
        self.lastStatus = step.getStatus()
        
        self.__updateStep(step)
        self._stepsDone.increment()
        self._store(self._stepsDone)
        
        self.info(magentaStr(step.getStatus().upper()) + ": %s, step %d" %
                  (step.funcName.get(), step._index))
        self.info("  %s" % step.endTime.datetime())
        if step.isFailed() and self.stepsExecutionMode == STEPS_PARALLEL:
            # In parallel mode the executor will exit to close
            # all working threads, so we need to close
            self._endRun()
        return doContinue
    
    def _stepsCheck(self):
        pass

    def _runSteps(self, startIndex):
        """ Run all steps defined in self._steps. """
        self._stepsDone.set(startIndex)
        self._numberOfSteps.set(len(self._steps))
        self.setRunning()
        self._originalRunMode = self.runMode.get()  # Keep the original value to set in sub-protocols
        self.runMode.set(MODE_RESUME)  # Always set to resume, even if set to restart
        self._store()
        
        if startIndex == len(self._steps):
            self.lastStatus = STATUS_FINISHED
            self.info("All steps seems to be FINISHED, nothing to be done.")
        else:
            self.lastStatus = self.status.get()
            self._stepsExecutor.runSteps(self._steps, 
                                         self._stepStarted, 
                                         self._stepFinished,
                                         self._stepsCheck)
        self.setStatus(self.lastStatus)
        self._store(self.status)
        
    def __deleteOutputs(self):
        """ This function should only be used from RESTART.
        It will remove output attributes from mapper and object.
        """
        attributes = [a[0] for a in self.iterOutputEM()]
        
        for attrName in attributes:
            attr = getattr(self, attrName)
            self.mapper.delete(attr)
            self.deleteAttribute(attrName)
            
        self._outputs.clear()
        self.mapper.store(self._outputs)
        
    def findAttributeName(self, attr):
        for attrName, attr in self.iterOutputEM():
            if attr.getObjId() == attr.getObjId():
                return attrName
        return None            
        
    def deleteOutput(self, output):
        attrName = self.findAttributeName(output)
        self.mapper.delete(output)
        self.deleteAttribute(attrName)
        if attrName in self._outputs:
            self._outputs.remove(attrName)
        self.mapper.store(self._outputs)
        self.mapper.commit()
        
    def __copyRelations(self, other):
        """ This will copy relations from protocol other to self """
        pass
    
    def copy(self, other, copyId=True):
        copyDict = Object.copy(self, other, copyId)
        self._store()
        self.mapper.deleteRelations(self)

        for r in other.getRelations():
            rName = r['name']
            rCreator = r['parent_id']
            rParent = r['object_parent_id']
            rChild = r['object_child_id']
            rParentExt = r['object_parent_extended']
            rChildExt = r['object_child_extended']

            if rParent in copyDict:
                rParent = copyDict.get(rParent).getObjId()

            if rChild in copyDict:
                rChild = copyDict.get(rChild).getObjId()

            self.mapper.insertRelationData(rName, rCreator, rParent, rChild,
                                           rParentExt, rChildExt)
        
    def getRelations(self):
        """ Return the relations created by this protocol. """
        return self.mapper.getRelationsByCreator(self)  
    
    def _defineRelation(self, relName, parentObj, childObj):
        """ Insert a new relation in the mapper using self as creator. """
        parentExt = None
        childExt = None
        
        if parentObj.isPointer():
            parentExt = parentObj.getExtended()
            parentObj = parentObj.getObjValue()
            
        if childObj.isPointer():
            childExt = childObj.getExtended()
            childObj = childObj.getObjValue()
            
        self.mapper.insertRelation(relName, self, parentObj, childObj,
                                   parentExt, childExt)
        
    def makePathsAndClean(self):
        """ Create the necessary path or clean
        if in RESTART mode. 
        """
        # Clean working path if in RESTART mode
        paths = [self._getPath(), self._getExtraPath(), self._getTmpPath(), self._getLogsPath()]
        
        if self.runMode == MODE_RESTART:
            cleanPath(*paths)
            self.__deleteOutputs()
            # Delete the relations created by this protocol
            # (delete this in both project and protocol db)
            self.mapper.deleteRelations(self)
        # Create workingDir, extra and tmp paths
        makePath(*paths)
        
    def cleanTmp(self):
        """ Delete all files and subdirectories under Tmp folder. """
        cleanPattern(self._getTmpPath('*'))
    
    def _run(self):
        # Check that a proper Steps executor have been set
        if self._stepsExecutor is None:
            raise Exception('Protocol.run: Steps executor should be set before running protocol')
        # Check the parameters are correct
        errors = self.validate()
        if len(errors):
            raise Exception('Protocol.run: Validation errors:\n' + '\n'.join(errors))
        
        #self.__backupSteps() # Prevent from overriden previous stored steps
        self._insertAllSteps() # Define steps for execute later
        #self._makePathsAndClean() This is done now in project
        startIndex = self.__findStartingStep() # Find at which step we need to start
        self.info(" Starting at step: %d" % (startIndex + 1))
        self._storeSteps() 
        self.info(" Running steps ")
        self._runSteps(startIndex)
    
    def _getEnviron(self):
        """ This function should return an environ variable
        that will be used when running new programs.
        By default, the protocol will use the one defined
        in the package that it belongs or None.
        """
        return getattr(self.getClassPackage(), '_environ', None)
        
    def runJob(self, program, arguments, **kwargs):
        if self.stepsExecutionMode == STEPS_SERIAL:
            kwargs['numberOfMpi'] = kwargs.get('numberOfMpi', self.numberOfMpi.get())
            kwargs['numberOfThreads'] = kwargs.get('numberOfThreads', self.numberOfThreads.get())
        else:
            kwargs['numberOfMpi'] = kwargs.get('numberOfMpi', 1)
            kwargs['numberOfThreads'] = kwargs.get('numberOfThreads', 1)
        if 'env' not in kwargs:
            #self._log.info("calling self._getEnviron...")
            kwargs['env'] = self._getEnviron()

        #self._log.info("runJob: cwd = %s" % kwargs.get('cwd', ''))
        #self._log.info("runJob: env = %s " % str(kwargs['env']))

        self._stepsExecutor.runJob(self._log, program, arguments, **kwargs)

    def run(self):
        """ Before calling this method, the working dir for the protocol
        to run should exists. 
        """
        self.__initLogs()
        
        self.info(greenStr('RUNNING PROTOCOL -----------------'))
        self._pid.set(os.getpid())
        self.info('          PID: %s' % self._pid)
        self.info('      Scipion: %s' % os.environ['SCIPION_VERSION'])
        self.info('   currentDir: %s' % os.getcwd())
        self.info('   workingDir: %s' % self.workingDir)
        self.info('      runMode: %s' % MODE_CHOICES[self.runMode.get()])
        try:
            self.info('          MPI: %d' % self.numberOfMpi)
            self.info('      threads: %d' % self.numberOfThreads)
        except Exception as e:
            self.info('  * Cannot get information about MPI/threads (%s)' % e)

        Step.run(self)
        self._endRun()

    def _endRun(self):
        """ Print some ending message and close some files. """   
        #self._store()
        self._store(self.summaryVar)
        self._store(self.methodsVar)
        self._store(self.endTime)

        if envVarOn('SCIPION_DEBUG_NOCLEAN'):
            self.warning('Not cleaning temporarly files since SCIPION_DEBUG_NOCLEAN is set to True.')
        elif not self.isFailed():
            self.info('Cleaning temporarly files....')
            self.cleanTmp()
            
        self.info(greenStr('------------------- PROTOCOL ' +
                           self.getStatusMessage().upper()))
        self.__closeLogs()
        
    def __initLogs(self):
        """ Open the log file overwriting its content if the protocol is going to be execute from zero. 
        Otherwise append the new content to the old one.
        Also open logs files and redirect the systems streams.
        """
        self._log = ScipionLogger(self.getLogPaths()[2]) 
               
        if self.runMode.get() == MODE_RESTART:
            mode = 'w+'
        else:     
            mode = 'a'
        # Backup the the system streams
        self.__stdErr = sys.stderr
        self.__stdOut = sys.stdout
        self.__openLogsFiles(mode)
        # Redirect the system streams to the protocol files
        sys.stdout = self.__fOut
        if envVarOn('SCIPION_SEPARATE_STDERR'):
            sys.stderr = self.__fErr
        else:
            sys.stderr = self.__fOut  # send stderr to wherever stdout appears
            # TODO: maybe do not show separate "run.stdout" and "run.stderr"
            # in the "Output Log" section if we put them together.

    def getLogPaths(self):
        return map(self._getLogsPath, ['run.stdout', 'run.stderr', 'run.log'])

    def getSteps(self):
        """ Return the steps.sqlite file under logs directory. """
        return self._steps
    
    def getStepsFile(self):
        """ Return the steps.sqlite file under logs directory. """
        return self._getLogsPath('steps.sqlite')

    def __openLogsFiles(self, mode):
        self.__fOut = open(self.getLogPaths()[0], mode)
        self.__fErr = open(self.getLogPaths()[1], mode)

    def __closeLogsFiles(self):
        self.__fOut.close()
        self.__fErr.close()
    
    def __closeLogs(self):
        self._log.close()
        # Restore system streams
        sys.stderr = self.__stdErr
        sys.stdout = self.__stdOut
        self.__closeLogsFiles()
        
    def _addChunk(self, txt, fmt=None):
        """
        Add text txt to self._buffer, with format fmt.
        fmt can be a color (like 'red') or a link that looks like 'link:url'.
        """
        # Make the text html-safe first.
        for x,y in [('&', 'amp'), ('<', 'lt'), ('>', 'gt')]:
            txt = txt.replace(x, '&%s;' % y)

        if fmt is None:
            self._buffer += txt
        elif fmt.startswith('link:'):
            url = fmt[len('link:'):]
            # Add the url in the TWiki style
            if url.startswith('http://'):
                self._buffer += '[[%s][%s]]' % (url, txt)
            else:
                from pyworkflow.web.pages import settings as django_settings
                absolute_url = django_settings.ABSOLUTE_URL
                self._buffer += '[[%s/get_log/?path=%s][%s]]' % (absolute_url, url, txt)
        else:
            self._buffer += '<font color="%s">%s</font>' % (fmt, txt)

    def getLogsAsStrings(self):

        outputs = []
        for fname in self.getLogPaths():
            if exists(fname):
                self._buffer = ''
                renderTextFile(fname, self._addChunk)
                outputs.append(self._buffer)
            else:
                outputs.append('File "%s" does not exist' % fname)
        return outputs

    def warning(self, message, redirectStandard=True):
        self._log.warning(message, redirectStandard)
        
    def info(self, message, redirectStandard=True):
        self._log.info(message, redirectStandard)
        
    def error(self, message, redirectStandard=True):
        self._log.error(message, redirectStandard)

    def debug(self, message):
        if envVarOn('SCIPION_DEBUG'):
            self.info(message)

    def getWorkingDir(self):
        return self.workingDir.get()
    
    def setWorkingDir(self, path):
        self.workingDir.set(path)

    def setMapper(self, mapper):
        """ Set a new mapper for the protocol to persist state. """
        self.mapper = mapper
        
    def getMapper(self):
        return self.mapper
        
    def getDbPath(self):
        return self._getLogsPath('run.db')
            
    def setStepsExecutor(self, executor=None):
        if executor is None:
            executor = StepExecutor(self.getHostConfig())
            
        self._stepsExecutor = executor
                
    def getFiles(self):
        resultFiles = set()
        for paramName, _ in self.getDefinition().iterPointerParams():
            attrPointer = getattr(self, paramName) # Get all self attribute that are pointers
            obj = attrPointer.get() # Get object pointer by the attribute
            if hasattr(obj, 'getFiles'):
                resultFiles.update(obj.getFiles()) # Add files if any
        return resultFiles | getFiles(self.workingDir.get())

    def getHostName(self):
        """ Get the execution host name """
        return self.hostName.get()
    
    def setHostName(self, hostName):
        """ Set the execution host name """ 
        self.hostName.set(hostName)
        
    def getHostConfig(self):
        """ Return the configuration host. """
        return self.hostConfig
    
    def setHostConfig(self, config):
        self.hostConfig = config
        # Never store the host config as part of the protocol, it is kept
        # in the configuration information, the hostname is enough
        self.hostConfig.setStore(False)
        
    def getJobId(self):
        """ Return the jobId associated to a running protocol. """
        return self._jobId.get()
    
    def setJobId(self, jobId):
        self._jobId.set(jobId)
        
    def getPid(self):
        return self._pid.get()
        
    def getRunName(self):
        runName = self.getObjLabel().strip()
        if not len(runName):
            runName = self.getDefaultRunName()
        return runName
    
    def getDefaultRunName(self):
        return '%s.%s' % (self.getClassName(), self.strId())
    
    @classmethod
    def getClassPackage(cls):
        """ Return the package module to which this protocol belongs
        """
        import pyworkflow.em as em
        em.getProtocols() # make sure the _package is set for each Protocol class
        return getattr(cls, '_package', scipion)
        
    @classmethod 
    def getClassPackageName(cls):
        return cls.getClassPackage().__name__.replace('pyworkflow.protocol.scipion', 'scipion')

    @classmethod
    def validatePackageVersion(cls, varName, errors):
        """
        Function to validate the the package version specified in configuration file
        ~/.config/scipion/scipion.conf is among the available options and it is
        properly installed.
        Params:
            package: the package object (ej: eman2 or relion). Package should contain the
                     following methods: getVersion(), getSupportedVersions()
            varName: the expected environment var containing the path (and version)
            errors: list to added error if found
        """
        package = cls.getClassPackage()
        packageName = cls.getClassPackageName()
        varValue = os.environ[varName]
        versions = ','.join(package.getSupportedVersions())

        errorMsg = None

        if not package.getVersion():
            errors.append("We could not detect *%s* version. " % packageName)
            errorMsg = "The path value should contains a valid version (%s)." % versions
        elif not os.path.exists(varValue):
            errors.append("Path of %s does not exists." % varName)
            errorMsg = "Check installed packages and versions with command:\n *scipion install --help*"

        if errorMsg:
            errors.append("%s = %s" % (varName, varValue))
            errors.append("Please, modify %s value in the configuration file:" % varName)
            errors.append("*~/.config/scipion/scipion.conf*")
            errors.append(errorMsg)
            errors.append("After fixed, you NEED TO RESTART THE PROJECT WINDOW")
        
    @classmethod
    def getClassLabel(cls):
        """ Return a more readable string representing the protocol class """
        label = cls.__dict__.get('_label', cls.__name__)
        label = "%s - %s" % (cls.getClassPackageName(), label)
        return label
    
    @classmethod
    def isBase(cls):
        """ Return True if this Protocol is a base class.
        Base classes should be marked with _label = None.
        """
        return cls.__dict__.get('_label', None) is None
        
    def getSubmitDict(self):
        """ Return a dictionary with the necessary keys to
        launch the job to a queue system.
        """
        queueName, queueParams = self.getQueueParams()
        hc = self.getHostConfig()
        
        script = self._getLogsPath(hc.getSubmitPrefix() + self.strId() + '.job')
        d = {'JOB_SCRIPT': script,
             'JOB_NODEFILE': script.replace('.job', '.nodefile'),
             'JOB_NAME': self.strId(),
             'JOB_QUEUE': queueName,
             'JOB_NODES': self.numberOfMpi.get(),
             'JOB_THREADS': self.numberOfThreads.get(),
             'JOB_CORES': self.numberOfMpi.get() * self.numberOfThreads.get(),
             'JOB_HOURS': 72,
             }
        d.update(queueParams)
        return d
    
    def useQueue(self):
        """ Return True if the protocol should be launched through a queue. """
        return self._useQueue.get()
    
    def getQueueParams(self):
        if self._queueParams.hasValue():
            return json.loads(self._queueParams.get())
        else:
            return '', {}
    
    def setQueueParams(self, queueParams):
        self._queueParams.set(json.dumps(queueParams))

    @property
    def numberOfSteps(self):
        return self._numberOfSteps.get(0)

    @property
    def stepsDone(self):
        """ Return the number of steps executed. """
        return self._stepsDone.get(0)

    def updateSteps(self):
        """ After the steps list is modified, this methods will update steps information.
        It will save the steps list and also the number of steps.
        """
        self._storeSteps()
        self._numberOfSteps.set(len(self._steps))
        self._store(self._numberOfSteps)

    def getStatusMessage(self):
        """ Return the status string and if running the steps done. 
        """
        msg = self.getStatus()
        if self.isRunning() or self.isAborted() or self.isFailed():
            msg += " (done %d/%d)" % (self.stepsDone, self.numberOfSteps)

        return msg
    
    def getRunMode(self):
        """ Return the mode of execution, either: MODE_RESTART or MODE_RESUME. """
        return self.runMode.get()
    
    # Methods that should be implemented in subclasses
    def _validate(self):
        """ This function can be overwritten by subclasses.
        Used from the public validate function.   
        """
        return []
     
    def validate(self):
        """ Check that input parameters are correct.
        Return a list with errors, if the list is empty, all was ok.         
        """
        errors = []
        # Validate that all input pointer parameters have a value
        for paramName, param in self.getDefinition().iterParams():
            attr = getattr(self, paramName) # Get all self attribute that are pointers
            paramErrors = []
            condition = self.evalParamCondition(paramName)
            if attr.isPointer():
                obj = attr.get()
                if condition and obj is None and not param.allowsNull:
                    paramErrors.append('cannot be EMPTY.')
            else:
                if condition:
                    paramErrors = param.validate(attr.get())
            label = param.label.get()
            errors += ['*%s* %s' % (label, err) for err in paramErrors]                
        # Validate specific for the subclass 
        childErrors = self._validate()        
        if childErrors:
            errors += childErrors
        
        return errors 
    
    def _warnings(self):
        """ Should be implemented in subclasses. See warning. """
        return []
    
    def warnings(self):
        """ Return some message warnings that can be errors.
        User should approve to execute a protocol with warnings. """
        return self._warnings()
    
    def _summary(self):
        """ Should be implemented in subclasses. See summary. """
        return ["No summary information."]
        
    def summary(self):
        """ Return a summary message to provide some information to users. """
        try:
            baseSummary = self._summary() or ['No summary information.']
        except Exception as ex:
            baseSummary = [str(ex)]        
            
        comments = self.getObjComment()
        if comments:
            baseSummary += ['', '*COMMENTS:* ', comments]
            
        if self.getError().hasValue():
            baseSummary += ['', '*ERROR:*', self.getError().get()]
            
        return baseSummary
    
    def getFileTag(self, fn):
        return "[[%s]]" % fn
    
    def getObjectTag(self, objName):
        if isinstance(objName, basestring):
            obj = getattr(self, objName, None)
        else:
            obj = objName
        
        if obj is None:
            return '*None*'
        
        if obj.isPointer():
            obj = obj.get() # get the pointed object
            if obj is None:
                return '*None*'
        
        return "[[sci-open:%s][%s]]" % (obj.getObjId(), obj.getNameId())
    #    return "[[javascript:launchViewer(%s)][%s]]" % (obj.getObjId(), obj.getNameId())
    
    def _citations(self):
        """ Should be implemented in subclasses. See citations. """
        return getattr(self, "_references", [])
    
    def __getPackageBibTex(self):
        """ Return the _bibtex from the package . """
        return getattr(self.getClassPackage(), "_bibtex", {})
    
    def _getCite(self, citeStr):
        bibtex = self.__getPackageBibTex()
        if citeStr in bibtex:
            text = self._getCiteText(bibtex[citeStr])
        else:
            text = "Reference with key *%s* not found." % citeStr
        return text
    
    def _getCiteText(self, cite, useKeyLabel=False):
        try:
            # Get the first author surname
            if useKeyLabel:
                label = cite['id']
            else:
                label = cite['author'].split(' and ')[0].split(',')[0].strip()
                label += ', et.al, %s, %s' % (cite['journal'], cite['year'])
            
            return '[[%s][%s]] ' % (cite['doi'].strip(), label)
        except Exception, ex:
            print "Error with citation: " + label
            print ex
            text = "Error with citation *%s*." % label
        return text
    
    def __getCitations(self, citations):
        """ From the list of citations keys, obtains the full
        info from the package _bibtex dict. 
        """
        bibtex = self.__getPackageBibTex()
        newCitations = []
        for c in citations:
            if c in bibtex:
                newCitations.append(self._getCiteText(bibtex[c]))
            else:
                newCitations.append(c)
        return newCitations
    
    def __getCitationsDict(self, citationList):
        """ Return a dictionary with Cite keys and the citation links. """
        bibtex = self.__getPackageBibTex()
        od = OrderedDict()
        for c in citationList:
            if c in bibtex:
                od[c] = self._getCiteText(bibtex[c])
            else:
                od[c] = c
            
        return od
        
    def getCitations(self):
        return self.__getCitationsDict(self._citations() or [])
            
    def getPackageCitations(self):
        return self.__getCitationsDict(getattr(self.getClassPackage(), "_references", []))
    
    def citations(self):
        """ Return a citation message to provide some information to users. """
        citations = self.getCitations().values()
        if citations:
            citations.insert(0, '*References:* ')
            
        packageCitations = self.getPackageCitations().values()
        if packageCitations:
            citations.append('*Package References:*')
            citations += packageCitations   
        if not citations:
            return ['No references provided']
        return citations

    def _methods(self):
        """ Should be implemented in subclasses. See methods. """
        return ["No methods information."]
        
    def getParsedMethods(self):
        """ Get the _methods results and parse possible cites. """
        try:
            baseMethods = self._methods() or []
            bibtex = self.__getPackageBibTex()
            parsedMethods = []
            for m in baseMethods:
                for bibId, cite in bibtex.iteritems():
                    k = '[%s]' % bibId
                    link = self._getCiteText(cite, useKeyLabel=True)
                    m = m.replace(k, link)
                parsedMethods.append(m)
        except Exception, ex:
            parsedMethods = ['ERROR generating methods info: %s' % ex]
        
        return parsedMethods
        
    def methods(self):
        """ Return a description about methods about current protocol execution. """
        # TODO: Maybe store the methods and not computing all times??
        return self.getParsedMethods() + [''] + self.citations()
        
    def runProtocol(self, protocol):
        """ Setup another protocol to be run from a workflow. """
        name = protocol.getClassName() + protocol.strId()
        #protocol.setName(name)
        protocol.setWorkingDir(self._getPath(name))
        protocol.setMapper(self.mapper)
        self.hostConfig.setStore(False)
        protocol.setHostConfig(self.getHostConfig())
        protocol.runMode.set(self._originalRunMode)
        protocol.makePathsAndClean()
        protocol.setStepsExecutor(self._stepsExecutor)
        protocol.run()
        self._store() #TODO: check if this is needed
        
    def isChild(self):
        """ Return true if this protocol was invoked from a workflow (another protocol)"""
        return self.hasObjParentId()
    
    def getStepsGraph(self, refresh=True):
        """ Build a graph taking into account the dependencies between steps. """
        from pyworkflow.utils.graph import Graph
        g = Graph(rootName='PROTOCOL')
        root = g.getRoot()
        root.label = 'Protocol'
        stepsDict = {}
        steps = self.loadSteps()
        
        for i, s in enumerate(steps):
            index = s._index or (i+1)
            sid = str(index)
            n = g.createNode(sid)
            n.label = sid
            stepsDict[sid] = n
            if s._prerequisites.isEmpty():
                root.addChild(n)
            else:
                for p in s._prerequisites:
                    stepsDict[p].addChild(n)
        return g
    
    def closeMappers(self):
        """ Close the mappers of all output Sets. """
        for _, attr in self.iterOutputAttributes(Set):
            attr.close()
            
    def loadMappers(self):
        """ Open mapper connections from previous closed outputs. """
        for _, attr in self.iterOutputAttributes(Set):
            attr.load()

    def allowsDelete(self, obj):
        return False
    
    def legacyCheck(self):
        """ Hook defined to run some compatibility checks
        before display the protocol.
        """
        pass


class LegacyProtocol(Protocol):
    """ Special subclass of Protocol to be used when a protocol class
    is not found. It means that have been removed or it is in another
    development branch. In such, we will use the LegacyProtocol to
    simply store the parameters and inputs/outputs."""

    def __str__(self):
        return self.getObjLabel()
                

#---------- Helper functions related to Protocols --------------------

def runProtocolMain(projectPath, protDbPath, protId):
    """ Main entry point when a protocol will be executed.
    This function should be called when:
    scipion runprotocol ...
    Params:
        projectPath: the absolute path to the project directory.
        protDbPath: path to protocol db relative to projectPath
        protId: id of the protocol object in db.
    """
    # Enter to the project directory and load protocol from db
    protocol = getProtocolFromDb(projectPath, protDbPath, protId, chdir=True)
    
    hostConfig = protocol.getHostConfig()  
  
    # Create the steps executor
    executor = None
    if protocol.stepsExecutionMode == STEPS_PARALLEL:
        if protocol.numberOfMpi > 1:
            # Handle special case to execute in parallel
            # We run "scipion run pyworkflow/...mpirun.py blah" instead of
            # calling directly "$SCIPION_PYTHON ...mpirun.py blah", so that
            # when it runs on a MPI node, it *always* has the scipion env.
            params = ['runprotocol', 'pw_protocol_mpirun.py',
                      projectPath, protDbPath, protId]
            retcode = runJob(None, pw.getScipionScript(), params,
                             numberOfMpi=protocol.numberOfMpi.get(),
                             hostConfig=hostConfig)
            sys.exit(retcode)
        elif protocol.numberOfThreads > 1:
            executor = ThreadStepExecutor(hostConfig,
                                          protocol.numberOfThreads.get()-1) 
    if executor is None:
        executor = StepExecutor(hostConfig)
    protocol.setStepsExecutor(executor)
    # Finally run the protocol
    protocol.run()
    
    
def runProtocolMainMPI(projectPath, protDbPath, protId, mpiComm):
    """ This function only should be called after enter in runProtocolMain
    and the proper MPI scripts have been started...so no validations 
    will be made.
    """ 
    protocol = getProtocolFromDb(projectPath, protDbPath, protId, chdir=True)
    hostConfig = protocol.getHostConfig()
    # Create the steps executor
    executor = MPIStepExecutor(hostConfig, protocol.numberOfMpi.get()-1, mpiComm)
    
    protocol.setStepsExecutor(executor)
    # Finally run the protocol
    protocol.run()        
    
     
def getProtocolFromDb(projectPath, protDbPath, protId, chdir=False):
    """ Retrieve the Protocol object from a given .sqlite file
    and the protocol id.
    """
    # We need this import here because from Project is imported
    # all from protocol indirectly, so if move this to the top
    # we get an import error
    if not exists(projectPath):
        raise Exception("ERROR: project path '%s' does not exist. " % projectPath)
        sys.exit(1)
    
    fullDbPath = os.path.join(projectPath, protDbPath)
    
    if not exists(fullDbPath):
        raise Exception("ERROR: protocol database '%s' does not exist. " % fullDbPath)
        sys.exit(1)
        
    from pyworkflow.project import Project
    project = Project(projectPath)
    project.load(dbPath=os.path.join(projectPath, protDbPath), chdir=chdir,
                 loadAllConfig=False)     
    protocol = project.getProtocol(protId)
    return protocol

