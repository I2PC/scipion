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
import datetime as dt
import pickle
import time

from pyworkflow.object import OrderedObject, String, List, Integer, Boolean, CsvList
from pyworkflow.utils.path import replaceExt, makePath, join, existsPath, cleanPath, getFolderFiles
from pyworkflow.utils.log import *

STATUS_SAVED = "saved" # Parameters saved for later use
STATUS_LAUNCHED = "launched"  # launched to queue system, only usefull for protocols
STATUS_RUNNING = "running"    # currently executing
STATUS_FAILED = "failed"      # it have been failed
STATUS_FINISHED = "finished"  # successfully finished
STATUS_WAITING_APPROVAL = "waiting approval"    # waiting for user interaction

# Following steps are only used for steps parallel execution
STATUS_READY = "ready" # The step is ready for execution, i.e. all requirements are done
STATUS_WAITING_OTHERS = "waiting_others" # There are some prerequisites steps that are not done yet


class Step(OrderedObject):
    """ Basic execution unit.
    It should defines its Input, Output
    and define a run method.
    """
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
        self._inputs = []
        self._outputs = []
        self._prerequisites = CsvList() # which steps needs to be done first
        self.status = String()
        self.initTime = String()
        self.endTime = String()
        self.error = String()
        self.isInteractive = Boolean(False)
        
    def _storeAttributes(self, attrList, attrDict):
        """ Store all attributes in attrDict as 
        attributes of self, also store the key in attrList.
        """
        for key, value in attrDict.iteritems():
            attrList.append(key)
            setattr(self, key, value)
        
    def _defineInputs(self, **args):
        """ This function should be used to define
        those attributes considered as Input""" 
        self._storeAttributes(self._inputs, args)
        
    def _defineOutputs(self, **args):
        """ This function should be used to specify
        expected outputs""" 
        self._storeAttributes(self._outputs, args)
           
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
        
    def setFailed(self, msg):
        """ Set the run failed and store an error message. """
        self.status.set(STATUS_FAILED)
        self.error.set(msg)
        
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
    def __init__(self, funcName=None, *funcArgs, **args):
        """ Receive the function to execute and the 
        parameters to call it""" 
        Step.__init__(self)
        self.func = None # Function should be set before run
        self.funcName = String(funcName)
        self.funcArgs = funcArgs
        self.argsStr = String(pickle.dumps(funcArgs))
        self.isInteractive.set(args.get('isInteractive', False))
        
    def _run(self):
        resultFiles = self.func(*self.funcArgs)
        if resultFiles and len(resultFiles):
            missingFiles = existsPath(*resultFiles)
            if len(missingFiles):
                raise Exception('Missing filePaths: ' + ' '.join(missingFiles))
            self.resultFiles = String(pickle.dumps(resultFiles))
    
    def _postconditions(self):
        """ This type of Step, will simply check
        as postconditions that the result filePaths exists""" 
        if not hasattr(self, 'resultFiles'):
            return True
        filePaths = pickle.loads(self.resultFiles.get())

        return len(existsPath(*filePaths)) == 0
    
    def __eq__(self, other):
        """ Compare with other FunctionStep""" 
        
        print 'self.funcName', self.funcName, 'other.funcName', other.funcName 
        print 'self.funcArgs == other.funcArgs', self.funcArgs == other.funcArgs 
        print 'self.argsStr == other.argsStr', self.argsStr == other.argsStr
        return (self.funcName == other.funcName and
                #self.funcArgs == other.funcArgs and
                self.argsStr == other.argsStr)
        
    def __ne__(self, other):
        return not self.__eq__(other)
        
            
class RunJobStep(FunctionStep):
    """ This Step will wrapper the commonly used function runJob
    for launching specific programs with some parameters""" 
    def __init__(self, programName=None, arguments=None, resultFiles=[], **args):
        FunctionStep.__init__(self, 'self.runJob', programName, arguments)
        # Define the function that will do the job and return result filePaths
        self.func = self._runJob
        self.mpi = 1
        self.threads = 1
        
    def _runJob(self, programName, arguments):
        """ Wrap around runJob function""" 
        self.runJob(None, programName, arguments, 
               numberOfMpi=self.mpi, numberOfThreads=self.threads)
        #TODO: Add the option to return resultFiles
             

MODE_RESUME = 0
MODE_RESTART = 1
MODE_CONTINUE = 2

STEPS_SERIAL = 0
STEPS_PARALLEL = 1
         
LEVEL_NORMAL = 0
                
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
        self._createVarsFromDefinition(**args)
        # For non-parallel protocols mpi=1 and threads=1
        if not hasattr(self, 'numberOfMpi'):
            self.numberOfMpi = Integer(1)
        if not hasattr(self, 'numberOfThreads'):
            self.numberOfThreads = Integer(1)
        # Maybe this property can be inferred from the 
        # prerequisites of steps, but is easier to keep it
        self.stepsExecutionMode = STEPS_SERIAL
        # Host name
        self.hostName = None
        # Expert level
        self.expertLevel = Integer(args.get('expertLevel', LEVEL_NORMAL))
        
    def getDefinition(self):
        """ Access the protocol definition. """
        return self._definition
    
    def getDefinitionParam(self, paramName):
        """ Return a _definition param give its name. """
        return self._definition.getParam(paramName)
    
    def evalCondition(self, paramName):
        """ Eval if the condition of paramName in _definition
        is satified with the current values of the protocol attributes. 
        """
        return self._definition.evalCondition(self, paramName)
    
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
                var = param.paramClass(args.get(paramName, param.default.get()))
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
        self.mapper.insertChild(self, key, child)
        
    def _defineSteps(self):
        """ Define all the steps that will be executed. """
        pass
    
    def __insertStep(self, step, **args):
        """ Insert a new step in the list. 
        Posible values for **args:
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
        
    def _insertFunctionStep(self, funcName, *funcArgs, **args):
        """ Input params:
        funcName: the string name of the function to be run in the Step.
        *funcArgs: the variable list of arguments to pass to the function.
        **args: see __insertStep
        """
        step = FunctionStep(funcName, *funcArgs, **args)
        step.func = getattr(self, funcName)
        return self.__insertStep(step, **args)
        
    def _insertRunJobStep(self, progName, progArguments, resultFiles=[], **args):
        """ Insert an Step that will simple call runJob function
        **args: see __insertStep
        """
        step = RunJobStep(progName, progArguments, resultFiles)
        step.runJob = self.runJob
        if self.stepsExecutionMode == STEPS_SERIAL:
            step.mpi = self.numberOfMpi.get()
            step.threads = self.numberOfThreads.get()
        return self.__insertStep(step, **args)
        
    def _enterWorkingDir(self):
        """ Change to the protocol working dir. """
        os.chdir(self.workingDir.get())
        
    def _leaveWorkingDir(self):
        """ This funcion make sense to use in conjunction 
        with _enterWorkingDir to go back to execution path.
        """
        os.chdir(self._currentDir)
        
    def __backupSteps(self):
        """ Store the Steps list in another variable to prevent
        overriden of stored steps when calling _defineSteps function.
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
        self._log.info("len(steps) " + str(len(self._steps)) + " len(prevSteps) " + str(len(self._prevSteps)))
        
        for i in range(n):
            newStep = self._steps[i]
            oldStep = self._prevSteps[i]
            print "i: ", i
            print " oldStep.status: ", str(oldStep.status)
            print " oldStep!=newStep", oldStep != newStep
            print " not post: ", not oldStep._postconditions()
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
        self._log.info("STARTED: " + step.funcName.get())
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
            self.setFailed("Protocol failed: " + step.error.get())
            self._log.error("Protocol failed: " + step.error.get())
        self.lastStatus = step.status.get()
        self._log.info("FINISHED: " + step.funcName.get())
        return doContinue
    
    def _runSteps(self, startIndex):
        """ Run all steps defined in self._steps. """
        self._steps.setStore(True) # Set steps to be stored
        self.setRunning()
        self._store()
        
        self.lastStatus = self.status.get()
        #status = STATUS_FINISHED # Just for the case doesn't enter in the loop
        
        self._stepsExecutor.runSteps(self._steps[startIndex:], 
                                     self._stepStarted, self._stepFinished)
        self.status.set(self.lastStatus)
        self._store(self.status)
        
    def _makePathsAndClean(self):
        """ Create the necessary path or clean
        if in RESTART mode. 
        """
        # Clean working path if in RESTART mode
        paths = [self.workingDir.get(), self._getExtraPath(), self._getTmpPath()]
        if self.runMode == MODE_RESTART:
            cleanPath(*paths)
            self.mapper.deleteChilds(self) # Clean old values
            self.mapper.insertChilds(self)
        # Create workingDir, extra and tmp paths
        makePath(*paths)        
    
    def _run(self):
        self.runJob = self._stepsExecutor.runJob
        self.__backupSteps() # Prevent from overriden previous stored steps
        self._defineSteps() # Define steps for execute later
        self._makePathsAndClean()
        startIndex = self.__findStartingStep() # Find at which step we need to start
        self.__cleanStepsFrom(startIndex) # 
        self._runSteps(startIndex)
        
    def run(self):
        self._log = self.__getLogger()
        self._log.info('RUNNING PROTOCOL -----------------')
        self._log.info('      runMode: %d' % self.runMode.get())
        self._log.info('   workingDir: ' + self.workingDir.get())
        #self.namePrefix = replaceExt(self._steps.getName(), self._steps.strId()) #keep
        self._currentDir = os.getcwd() 
        self._run()
        outputs = [getattr(self, o) for o in self._outputs]
        #self._store(self.status, self.initTime, self.endTime, *outputs)
        self._store()
        self._log.info('------------------- PROTOCOL FINISHED')
        
    def __getLogger(self):
        #Initialize log
        logFile = self._getPath('log', 'protocol.log')
        return getFileLogger(logFile)

    def getFiles(self):
        resultFiles = set()
        for paramName, _ in self.getDefinition().iterPointerParams():
            attrPointer = getattr(self, paramName) # Get all self attribute that are pointers
            obj = attrPointer.get() # Get object pointer by the attribute
            if hasattr(obj, 'getFiles'):
                resultFiles.update(obj.getFiles()) # Add files if any
        return resultFiles | getFolderFiles(self.workingDir.get())

    def getHostName(self):
        """ Get the execution host name """
        return self.hostName
    
    def setHostName(self, hostName):
        """ Set the execution host name """ 
        self.hostName = hostName
        
    def getElapsedTime(self):
        """ Return the time that protocols
        took to run (or the actual running time 
        if still is running )
        """
        f = "%Y-%m-%d %H:%M:%S.%f"
        t1 = dt.datetime.strptime(self.initTime.get(), f)
        if self.status == STATUS_RUNNING:
            t2 = dt.datetime.now()
        else:
            t2 = dt.datetime.strptime(self.endTime.get(), f)
        elapsed = t2 - t1
        
        return elapsed
    
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
        return self._validate()
    
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
        return self._summary()
        
