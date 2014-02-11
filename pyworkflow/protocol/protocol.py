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
from pyworkflow.protocol.constants import STATUS_RUNNING, STATUS_FAILED
"""
This modules contains classes required for the workflow
execution and tracking like: Step and Protocol
"""

import os
import sys
import datetime as dt
import pickle
import time

from pyworkflow.object import *
from pyworkflow.utils.path import replaceExt, makeFilePath, join, missingPaths, cleanPath, getFiles
from pyworkflow.utils.log import *
from executor import StepExecutor, ThreadStepExecutor, MPIStepExecutor
from constants import *
from params import Form



class Step(OrderedObject):
    """ Basic execution unit.
    It should defines its Input, Output
    and define a run method.
    """
    
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
        self._prerequisites = CsvList() # which steps needs to be done first
        self.status = String()
        self.initTime = String()
        self.endTime = String()
        self._error = String()
        self.isInteractive = Boolean(False)
           
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
        
    def setFailed(self, msg):
        """ Set the run failed and store an error message. """
        self._error.set(msg)
        self.status.set(STATUS_FAILED)
        
    def setAborted(self):
        """ Set the status to aborted and updated the endTime. """
        self.endTime.set(dt.datetime.now())
        self._error.set("Aborted by user.")
        self.status.set(STATUS_ABORTED)

    def getStatus(self):
        status = self.status.get()
        if status == None:
            status= "new"
        return status
    
    def setStatus(self, value):
        return self.status.set(value)
    
    def isActive(self):
        return self.getStatus() in ACTIVE_STATUS
    
    def isFinished(self):
        return self.getStatus() == STATUS_FINISHED
    
    def isRunning(self):
        return self.getStatus() == STATUS_RUNNING
    
    def isFailed(self):
        return self.getStatus() == STATUS_FAILED
    
    def run(self):
        """ Do the job of this step"""
        self.setRunning() 
        try:
            self._run()
            if self.status == STATUS_RUNNING:
                if self.isInteractive.get():
                    # If the Step is interactive, after run
                    # it will be waiting for use to mark it as DONE
                    status = STATUS_WAITING_APPROVAL
                else:
                    status = STATUS_FINISHED
                self.status.set(status)
        except Exception, e:
            self.setFailed(str(e))
            import traceback
            traceback.print_exc()            
            #raise #only in development
        finally:
            self.endTime.set(dt.datetime.now())


class FunctionStep(Step):
    """ This is a Step wrapper around a normal function
    This class will ease the insertion of Protocol function steps
    throught the function _insertFunctionStep""" 
    
    def __init__(self, func=None, funcName=None, *funcArgs, **args):
        """ 
         Params:
            func: the function that will be executed.
            funcName: the name assigned to that function (will be stored)
            *funcArgs: argument list passed to the function (serialized and stored)
            **args: extra parameters.
        """ 
        Step.__init__(self)
        self._func = func # Function should be set before run
        self._args = funcArgs
        self.funcName = String(funcName)
        self.argsStr = String(pickle.dumps(funcArgs))
        self.isInteractive.set(args.get('isInteractive', False))
        
    def _runFunc(self):
        """ Return the possible result files after running the function. """
        return self._func(*self._args)
    
    def _run(self):
        """ Run the function and check the result files if any. """
        resultFiles = self._runFunc() 
        if resultFiles and len(resultFiles):
            missingFiles = missingPaths(*resultFiles)
            if len(missingFiles):
                raise Exception('Missing filePaths: ' + ' '.join(missingFiles))
            self.resultFiles = String(pickle.dumps(resultFiles))
    
    def _postconditions(self):
        """ This type of Step, will simply check
        as postconditions that the result filePaths exists""" 
        if not hasattr(self, 'resultFiles'):
            return True
        filePaths = pickle.loads(self.resultFiles.get())

        return len(missingPaths(*filePaths)) == 0
    
    def __eq__(self, other):
        """ Compare with other FunctionStep""" 
        return (self.funcName == other.funcName and
                self.argsStr == other.argsStr)
        
    def __ne__(self, other):
        return not self.__eq__(other)
        
            
class RunJobStep(FunctionStep):
    """ This Step will wrapper the commonly used function runJob
    for launching specific programs with some parameters.
    The runJob function should be provided by the protocol
    when inserting a new RunJobStep"""
     
    def __init__(self, runJobFunc=None, programName=None, arguments=None, resultFiles=[], **args):
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
             
                
class Protocol(Step):
    """ The Protocol is a higher type of Step.
    It also have the inputs, outputs and other Steps properties,
    but contains a list of steps that are executed
    """
    
    def __init__(self, **args):
        Step.__init__(self, **args)        
        self._steps = List() # List of steps that will be executed
        self.workingDir = String(args.get('workingDir', '.')) # All generated filePaths should be inside workingDir
        self.mapper = args.get('mapper', None)
        self._inputs = []
        self._outputs = CsvList()
        self._definition = Form()
        self._defineParams(self._definition)
        self._createVarsFromDefinition(**args)
        self.__stdOut = None
        self.__stdErr = None
        self.__fOut = None
        self.__fErr = None
        
        # For non-parallel protocols mpi=1 and threads=1
        if not hasattr(self, 'numberOfMpi'):
            self.numberOfMpi = Integer(1)
        
        if not hasattr(self, 'numberOfThreads'):
            self.numberOfThreads = Integer(1)
        
        # Check if MPI or threads are passed in **args, mainly used in tests
        if 'numberOfMpi' in args:
            self.numberOfMpi.set(args.get('numberOfMpi'))
        
        if 'numberOfThreads' in args:
            self.numberOfThreads.set(args.get('numberOfThreads'))            
        
        if not hasattr(self, 'hostName'):
            self.hostName = String(args.get('hostName', 'localhost'))
        
        # Maybe this property can be inferred from the 
        # prerequisites of steps, but is easier to keep it
        self.stepsExecutionMode = STEPS_SERIAL
        # Expert level
        self.expertLevel = Integer(args.get('expertLevel', LEVEL_NORMAL))
        # Run mode
        self.runMode = Integer(args.get('runMode', MODE_RESUME))
        # Use queue system?
        self._useQueue = Boolean()
        
        self._jobId = String() # Store queue job id
        self._pid = Integer()
        self._stepsExecutor = None
        
        # For visualization
        self.allowHeader = Boolean(True)        
        
    def _storeAttributes(self, attrList, attrDict):
        """ Store all attributes in attrDict as 
        attributes of self, also store the key in attrList.
        """
        for key, value in attrDict.iteritems():
            if key not in attrList:
                attrList.append(key)
            setattr(self, key, value)
        
    def _defineInputs(self, **args):
        """ This function should be used to define
        those attributes considered as Input""" 
        self._storeAttributes(self._inputs, args)
        
    def _defineOutputs(self, **args):
        """ This function should be used to specify
        expected outputs""" 
        
        for k, v in args.iteritems():
            
            if hasattr(self, k):
                self._deleteChild(k, v)
            self._insertChild(k, v)
            
        self._storeAttributes(self._outputs, args)
        
        
                
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
    
    def getDefinitionParam(self, paramName):
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
        return self.getDefinitionParam(paramName).choices[index]
    
    def evalParamCondition(self, paramName):
        """ Eval if the condition of paramName in _definition
        is satified with the current values of the protocol attributes. 
        """
        return self._definition.evalParamCondition(self, paramName)
    
    def evalExpertLevel(self, paramName):
        """ Return True if the param has an expert level is less than 
        the one for the whole protocol. 
        """
        return self.getDefinitionParam(paramName).expertLevel.get() <= self.expertLevel.get()
    
    def iterDefinitionAttributes(self):
        """ Iterate over all the attributes from definition. """
        for paramName, _ in self._definition.iterParams():
            yield paramName, getattr(self, paramName)
            
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
            if isinstance(attr, PointerList) and attr.hasValue():
                for item in attr:
                    yield key, item
            elif attr.isPointer() and attr.hasValue():
                yield key, attr
                
    def iterOutputAttributes(self, outputClass):
        """ Iterate over the outputs produced by this protocol. """
        for key, attr in self.getAttributes():
            if isinstance(attr, outputClass):
                yield key, attr
            
    def copyDefinitionAttributes(self, other):
        """ Copy definition attributes to other protocol. """
        for paramName, _ in self.iterDefinitionAttributes():
            self.copyAttributes(other, paramName)
        
    def _createVarsFromDefinition(self, **args):
        """ This function will setup the protocol instance variables
        from the Protocol Class definition, taking into account
        the variable type and default values.
        """
        if hasattr(self, '_definition'):
            for paramName, param in self._definition.iterParams():
                # Create the var with value comming from args or from 
                # the default param definition
                var = param.paramClass(value=args.get(paramName, param.default.get()))
                setattr(self, paramName, var)
        else:
            print "FIXME: Protocol '%s' has not DEFINITION" % self.getClassName()
        
    
    def _getFilename(self, key, **args):
        return self._templateDict[key] % args
    
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
       
        setattr(self, key, child)
        if self.hasObjId():
            self.mapper.insertChild(self, key, child)
        
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
    
    def __insertStep(self, step, **args):
        """ Insert a new step in the list. 
        Params:
         **args:
            prerequisites: a list with the steps index that need to be done 
                           previous than the current one."""
        prerequisites = args.get('prerequisites', None)
        
        if prerequisites is None:
            if len(self._steps):
                step._prerequisites.append(len(self._steps)) # By default add the previous step as prerequisite
        else:
            for i in prerequisites:
                step._prerequisites.append(i)
                
        self._steps.append(step)
        # Return step number
        return len(self._steps)
        
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
    
    def _getBasePath(self, path):
        """ Take the basename of the path and get the path
        relative to working dir of the protocol. 
        """
        return self._getPath(os.path.basename(path))
        
    def _insertFunctionStep(self, funcName, *funcArgs, **args):
        """ 
         Params:
           funcName: the string name of the function to be run in the Step.
           *funcArgs: the variable list of arguments to pass to the function.
           **args: see __insertStep
        """
        # Get the function give its name
        func = getattr(self, funcName, None)
        # Ensure the protocol instance have it and is callable
        if not func:
            raise Exception("Protocol._insertFunctionStep: '%s' function is not member of the protocol" % funcName)
        if not callable(func):
            raise Exception("Protocol._insertFunctionStep: '%s' is not callable" % funcName)
        step = FunctionStep(func, funcName, *funcArgs, **args)
        
        return self.__insertStep(step, **args)
        
    def _insertRunJobStep(self, progName, progArguments, resultFiles=[], **args):
        """ Insert an Step that will simple call runJob function
        **args: see __insertStep
        """
        self._insertFunctionStep('runJob', progName, progArguments)
            
    def _enterDir(self, path):
        """ Enter into a new directory path and store the current path.
        The current path will be used in _leaveDir, but nested _enterDir
        are not allowed since self._currentDir is overriden.
        """
        self._currentDir = os.getcwd()
        os.chdir(path)
        
    def _leaveDir(self): 
        """ This method should be called after a call to _enterDir
        to return to the previous location. 
        """
        os.chdir(self._currentDir)        
                      
    def _enterWorkingDir(self):
        """ Change to the protocol working dir. """
        self._enterDir(self.workingDir.get())
        
    def _leaveWorkingDir(self):
        """ This funcion make sense to use in conjunction 
        with _enterWorkingDir to go back to execution path.
        """
        self._leaveDir()
        
    def __backupSteps(self):
        """ Store the Steps list in another variable to prevent
        overriden of stored steps when calling _insertAllSteps function.
        This is need to later find in which Step will start the run
        if the RESUME mode is used.
        """
        self._steps.setStore(False)
        self._prevSteps = self._steps
        self._steps = List() # create a new object for steps
        
    def __findStartingStep(self):
        """ From a previous run, compare self._steps and self._prevSteps
        to find which steps we need to start at, skipping sucessful done 
        and not changed steps. Steps that needs to be done, will be deleted
        from the previous run storage.
        """
        if self.runMode == MODE_RESTART:
            return 0
        
        n = min(len(self._steps), len(self._prevSteps))
        self.info("len(steps) " + str(len(self._steps)) + " len(prevSteps) " + str(len(self._prevSteps)))
        
        for i in range(n):
            newStep = self._steps[i]
            oldStep = self._prevSteps[i]
#            print "i: ", i
#            print " oldStep.status: ", str(oldStep.status)
#            print " oldStep!=newStep", oldStep != newStep
#            print " not post: ", not oldStep._postconditions()
            if (oldStep.status.get() != STATUS_FINISHED or
                newStep != oldStep or 
                not oldStep._postconditions()):
                return i
            
        return n
    
    def __cleanStepsFrom(self, index):
        """ Clean from the persistence all steps in self._prevSteps
        from that index. After this function self._steps is updated
        with steps from self._prevSteps (<index) and self._prevSteps is 
        deleted since is no longer needed.
        """
        self._steps[:index] = self._prevSteps[:index]
        
        #for oldStep in self._prevSteps[index:]:
        for oldStep in self._prevSteps: # This delete all and insert them later
            self.mapper.delete(oldStep)
            
        self._prevSteps = []
        
    def _stepStarted(self, step):
        """This function will be called whenever an step
        has started running.
        """
        self.info("STARTED: " + step.funcName.get())
        self.status.set(step.status)
        self._store(step)
    
    def _stepFinished(self, step):
        """This function will be called whenever an step
        has finished its run.
        """
        self.endTime.set(step.endTime.get())
        self._store(self.endTime, step)
        doContinue = True
        if step.status == STATUS_WAITING_APPROVAL:
            doContinue = False
        elif step.status == STATUS_FAILED:
            doContinue = False
            self.setFailed("Protocol failed: " + step.getError().get())
            self.error("Protocol failed: " + step.getError().get())
        self.lastStatus = step.status.get()
        self.info("FINISHED: " + step.funcName.get())
        return doContinue
    
    def _runSteps(self, startIndex):
        """ Run all steps defined in self._steps. """
        self._steps.setStore(True) # Set steps to be stored
        self.setRunning()
        self._originalRunMode = self.runMode.get() # Keep the original value to set in sub-protocols
        self.runMode.set(MODE_RESUME) # Always set to resume, even if set to restart
        self._store()
        
        self.lastStatus = self.status.get()
        #status = STATUS_FINISHED # Just for the case doesn't enter in the loop
        
        self._stepsExecutor.runSteps(self._steps[startIndex:], 
                                     self._stepStarted, self._stepFinished)
        self.status.set(self.lastStatus)
        self._store(self.status)
        
    def __deleteOutputs(self):
        """ This function should only be used from RESTART.
        It will remove output attributes from mapper and object.
        """
        for o in self._outputs:
            if hasattr(self, o):
                self.mapper.delete(getattr(self, o))
                self.deleteAttribute(o)
            else:
                print "DEBUG: error, %s is not found in protocol" % o
            
        self._outputs.clear()
        self.mapper.store(self._outputs)
        
    def __copyRelations(self, other):
        """ This will copy relations from protocol other to self """
        pass
    
    def copy(self, other, copyId=True):
        copyDict = Object.copy(self, other, copyId)
        self._store()
        for r in other.getRelations():
            rName = r['name']
            rCreator = r['parent_id']
            rParent = r['object_parent_id']
            rChild = r['object_child_id']
            
            if rParent in copyDict:
                rParent = copyDict.get(rParent).getObjId()
                            
            if rChild in copyDict:
                rChild = copyDict.get(rChild).getObjId()
            
            self.mapper.insertRelationData(rName, rCreator, rParent, rChild)
        
        
    def getRelations(self):
        """ Return the relations created by this protocol. """
        return self.mapper.getRelationsByCreator(self)  
    
    def _defineRelation(self, relName, parentObj, childObj):
        """ Insert a new relation in the mapper using self as creator. """
        self.mapper.insertRelation(relName, self, parentObj, childObj)
        
    def makePathsAndClean(self):
        """ Create the necessary path or clean
        if in RESTART mode. 
        """
        # Clean working path if in RESTART mode
        paths = [self._getPath(), self._getExtraPath(), self._getTmpPath(), self._getLogsPath()]
        
        if self.runMode == MODE_RESTART:
            cleanPath(*paths)
            self.__deleteOutputs()
        # Create workingDir, extra and tmp paths
        makePath(*paths)
    
    def _run(self):
        # Check that a proper Steps executor have been set
        if self._stepsExecutor is None:
            raise Exception('Protocol.run: Steps executor should be set before running protocol')
        # Check the parameters are correct
        errors = self.validate()
        if len(errors):
            raise Exception('Protocol.run: Validation errors:\n' + '\n'.join(errors))
        
        self.__backupSteps() # Prevent from overriden previous stored steps
        self._insertAllSteps() # Define steps for execute later
        #self._makePathsAndClean() This is done now in project
        startIndex = self.__findStartingStep() # Find at which step we need to start
        self.info(" Starting at index: %d" % startIndex)
        self.__cleanStepsFrom(startIndex) # 
        self._runSteps(startIndex)
    
    def runJob(self, program, arguments, **kwargs):
        if self.stepsExecutionMode == STEPS_SERIAL:
            kwargs['numberOfMpi'] = kwargs.get('numberOfMpi', self.numberOfMpi.get())
            kwargs['numberOfThreads'] = kwargs.get('numberOfThreads', self.numberOfThreads.get())
        else:
            kwargs['numberOfMpi'] = kwargs.get('numberOfMpi', 1)
            kwargs['numberOfThreads'] = kwargs.get('numberOfThreads', 1)           
            
        self._stepsExecutor.runJob(self._log, program, arguments, **kwargs)
        
    def run(self):
        """ Before calling this method, the working dir for the protocol
        to run should exists. 
        """
        self.__initLogs()
        
        self.info('RUNNING PROTOCOL -----------------')
#        self.info('RUNNING PROTOCOL info -----------------')
#        self.warning('RUNNING PROTOCOL warning-----------------')
#        self.error('RUNNING PROTOCOL error-----------------')
#        self.info('RUNNING PROTOCOL info red-----------------', True)
#        self.warning('RUNNING PROTOCOL warning red-----------------', True)
#        self.error('RUNNING PROTOCOL error red-----------------', True)
#        self.info('        jobId: %s' % self.getJobId())
#        self.info('          pid: %s' % os.getpid())
#        self.info('         ppid: %s' % os.getppid())
        self._pid.set(os.getpid())
        self.info('   currentDir: %s' % os.getcwd())
        self.info('   workingDir: ' + self.workingDir.get())
        self.info('      runMode: %d' % self.runMode.get())
        #self.namePrefix = replaceExt(self._steps.getName(), self._steps.strId()) #keep
        #self._run()
        Step.run(self)
        
        #outputs = [getattr(self, o) for o in self._outputs]
        self._store()
        self.info('------------------- PROTOCOL FINISHED')
        self.__closeLogs()
        
    def __initLogs(self):
        """ Open the log file overwriting its content if the protocol is going to be execute from zero. 
        Otherwise append the new content to the old one.
        Also open logs files and redirect the systems streams.
        """
        self._log = ScipionLogger(self.__getLogPaths()[2]) 
               
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
        sys.stderr = self.__fErr
    
    def __getLogPaths(self):
        return self._getLogsPath('run.stdout'), self._getLogsPath('run.stderr'), self._getLogsPath('run.log')
    
    def __openLogsFiles(self, mode):
        self.__fOut = open(self.__getLogPaths()[0], mode)
        self.__fErr = open(self.__getLogPaths()[1], mode)
        
    def __closeLogsFiles(self):       
        self.__fOut.close()
        self.__fErr.close()
    
    def __closeLogs(self):
        self._log.close()
        # Restore system streams
        sys.stderr = self.__stdErr
        sys.stdout = self.__stdOut
        self.__closeLogsFiles()
        
    def getLogsAsStrings(self):
        fOutString = fErrString = fScpnString = ''
        if os.path.exists(self.__getLogPaths()[0]) and os.path.exists(self.__getLogPaths()[1]) and os.path.exists(self.__getLogPaths()[2]):
            self.__openLogsFiles('r')
            fOutString = self.__fOut.read()
            fErrString = self.__fErr.read()
            self.__closeLogsFiles()
            
            fScpn = open(self.__getLogPaths()[2], 'r')
            fScpnString = fScpn.read()
            fScpn.close()
            
        return fOutString, fErrString, fScpnString
    
    def warning(self, message, redirectStandard=True):
        self._log.warning(message, redirectStandard)
        
    def info(self, message, redirectStandard=True):
        self._log.info(message, redirectStandard)
        
    def error(self, message, redirectStandard=True):
        self._log.error(message, redirectStandard)
        
    def getWorkingDir(self):
        return self.workingDir.get()
    
    def setWorkingDir(self, path):
        self.workingDir.set(path)

    def setMapper(self, mapper):
        """ Set a new mapper for the protocol to persist state. """
        self.mapper = mapper
        
    def getDbPath(self):
        return self._getLogsPath('run.db')
            
    def setStepsExecutor(self, executor):
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
    def getClassLabel(cls):
        """ Return a more readable string representing the protocol class """
        label = getattr(cls, '_label', cls.__name__)
        label = "%s - %s" % (cls._package.__name__, label)
        return label
        
    def getSubmitDict(self):
        """ Return a dictionary with the necessary keys to
        launch the job to a queue system.
        """
        script = self._getLogsPath(self.strId() + '.job')
        d = {'JOB_SCRIPT': script,
             'JOB_NODEFILE': script.replace('.job', '.nodefile'),
             'JOB_NAME': self.strId(),
             'JOB_QUEUE': 'default',
             'JOB_NODES': self.numberOfMpi.get(),
             'JOB_THREADS': self.numberOfThreads.get(),
             'JOB_CORES': self.numberOfMpi.get() * self.numberOfThreads.get(),
             'JOB_HOURS': 72,
             }
        return d
    
    def useQueue(self):
        """ Return True if the protocol should be launched throught a queue. """
        return self._useQueue.get()
        
    def getElapsedTime(self):
        """ Return the time that protocols
        took to run (or the actual running time 
        if still is running )
        """
        elapsed = None
        if self.initTime.hasValue():
            f = "%Y-%m-%d %H:%M:%S.%f"
            t1 = dt.datetime.strptime(self.initTime.get(), f)
            if self.endTime.hasValue():
                t2 = dt.datetime.strptime(self.endTime.get(), f)
            else:
                t2 = dt.datetime.now()
            elapsed = t2 - t1
        
        return elapsed
    
    def getStepsDone(self):
        """ Return the number of steps executed. """
        done = 0
        for s in self._steps:
            if s.isFinished():
                done += 1
        return done
            
    def getStatusMessage(self):
        """ Return the status string and if running the steps done. 
        """
        msg = self.getStatus()
        if self.isRunning():
            done = self.getStepsDone()
            msg += " (done %d/%d)" % (done, len(self._steps))
        
        return msg
    
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
                if condition and obj is None and not param.allowNull:
                    paramErrors.append('cannot be EMPTY.')
            else:
                if condition:
                    paramErrors = param.validate(attr.get())
            label = param.label.get()
            errors += ['*%s* %s' % (label, err) for err in paramErrors]                
        # Validate specific for the subclass
        errors += self._validate()
        
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
        error = ''
        if self.getError().hasValue():
            error = '*ERROR:*\n' + self.getError().get()
        baseSummary = self._summary()
        if not baseSummary:
            baseSummary = []
        return baseSummary + ['', '*Comments:* ', self.getObjComment(), error]
    
    def _citations(self):
        """ Should be implemented in subclasses. See citations. """
        return getattr(self, "_references", [])
    
    def __getPackageBibTex(self):
        """ Return the _bibtex from the package . """
        return getattr(self._package, "_bibtex", {})
     
    def __getCiteText(self, cite, useKeyLabel=False):
        # Get the first author surname
        if useKeyLabel:
            label = cite['id']
        else:
            label = cite['author'].split('and')[0].split(',')[0].strip()
            label += ', et.al, %s, %s' % (cite['journal'], cite['year'])
        
        return '[[%s][%s]] ' % (cite['doi'].strip(), label)
    
    def __getCitations(self, citations):
        """ From the list of citations keys, obtains the full
        info from the package _bibtex dict. 
        """
        bibtex = self.__getPackageBibTex()
        newCitations = []
        for c in citations:
            if c in bibtex:
                newCitations.append(self.__getCiteText(bibtex[c]))
            else:
                newCitations.append(c)
        return newCitations
    
    def citations(self):
        """ Return a citation message to provide some information to users. """
        citations = self._citations()
        if not citations:
            citationsFinal = []
        else:
            citationsFinal = ['*References:* '] + citations

        packageCitations = getattr(self._package, "_references", [])
        if packageCitations:
            citationsFinal = citationsFinal +['*Package References:*'] + packageCitations   
            
        return self.__getCitations(citationsFinal)    

    def _methods(self):
        """ Should be implemented in subclasses. See methods. """
        return ["No methods information."]
        
    def methods_old(self):
        """ Return a description about methods about current protocol execution. """
        baseMethods = self._methods()
        if not baseMethods:
            baseMethods = []
        return baseMethods + [''] + self.citations()
        
    def methods(self):
        """ Return a description about methods about current protocol execution. """
        # TODO: Maybe store the methods and not computing all times??
        baseMethods = self._methods()
        if not baseMethods:
            baseMethods = []
        else:
            bibtex = self.__getPackageBibTex()
            parsedMethods = []
            for m in baseMethods:
                for bibId, cite in bibtex.iteritems():
                    k = '[%s]' % bibId
                    link = self.__getCiteText(cite, useKeyLabel=True)
                    m = m.replace(k, link)
                parsedMethods.append(m)
            baseMethods = parsedMethods
            
        return baseMethods + [''] + self.citations()
        
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
        
    def isChild(self):
        """ Return true if this protocol was invoked from a workflow (another protocol)"""
        return self.hasObjParentId()
        
                
#---------- Helper functions related to Protocols --------------------

def getProtocolFromDb(dbPath, protId, protDict):
    from pyworkflow.mapper import SqliteMapper
    mapper = SqliteMapper(dbPath, protDict)
    protocol = mapper.selectById(protId)
    if protocol is None:
        raise Exception("Not protocol found with id: %d" % protId)
    protocol.setMapper(mapper)
    return protocol
    
def runProtocolFromDb(dbPath, protId, protDict, mpiComm=None):
    """ Retrieve a protocol stored in a db and run it. 
    Params:
     path: file path to read the mapper
     protId: the protocolo id in the mapper
     protDict: the dictionary with protocol classes.
     mpiComm: MPI connection object (only used when executing with MPI)
    """
    protocol = getProtocolFromDb(dbPath, protId, protDict)
    hostConfig = protocol.getHostConfig()
    # Create the steps executor
    executor = None
    if protocol.stepsExecutionMode == STEPS_PARALLEL:
        if protocol.numberOfMpi > 1:
            if mpiComm is None:
                raise Exception('Trying to create MPIStepExecutor and mpiComm is None')
            executor = MPIStepExecutor(hostConfig,
                                       protocol.numberOfMpi.get(), mpiComm)
        elif protocol.numberOfThreads > 1:
            executor = ThreadStepExecutor(hostConfig,
                                          protocol.numberOfThreads.get()) 
    if executor is None:
        executor = StepExecutor(hostConfig)
    protocol.setStepsExecutor(executor)
    # Finally run the protocol
    protocol.run()        

