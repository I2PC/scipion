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
from utils.path import replaceExt
"""
This modules contains classes required for the workflow
execution and tracking like: Step and Protocol
"""

import datetime as dt
import pickle

from pyworkflow.object import FakedObject, String, List, Integer

STATUS_LAUNCHED = "launched"  # launched to queue system
STATUS_RUNNING = "running"    # currently executing
STATUS_FAILED = "failed"      # it have been failed
STATUS_FINISHED = "finished"  # successfully finished
STATUS_WAITING = "waiting"    # waiting for user interaction

class Step(FakedObject):
    """Basic execution unit.
    It should defines its Input, Output
    and define a run method"""
    def __init__(self, **args):
        FakedObject.__init__(self, **args)
        self._inputs = []
        self._outputs = []
        self.addAttribute('status', String)
        self.addAttribute('initTime', String)
        self.addAttribute('endTime', String)
        self.addAttribute('error', String)
        
    def _storeAttributes(self, attrList, attrDict):
        """Store all attributes in attrDict as 
        attributes of self, also store the key in attrList"""
        for key, value in attrDict.iteritems():
            attrList.append(key)
            setattr(self, key, value)
        
    def defineInputs(self, **args):
        """This function should be used to define
        those attributes considered as Input"""
        self._storeAttributes(self._inputs, args)
        
    def defineOutputs(self, **args):
        """This function should be used to specify
        expected outputs"""
        self._storeAttributes(self._outputs, args)
    
    def preconditions(self):
        """Check if the necessary conditions to
        step execution are met"""
        return True
    
    def postconditions(self):
        """Check if the step have done well its task
        and accomplish its results"""
        return True
    
    def _run(self):
        """This is the function that will do the real job.
        It should be override by sub-classes."""
        pass
    
    def run(self):
        """Do the job of this step"""
        self.initTime = str(dt.datetime.now())
        self.endTime = None
        try:
            self._run()
            self.status = STATUS_FINISHED
        except Exception, e:
            self.status = STATUS_FAILED
            self.error = str(e)
            raise #only in development
        finally:
            self.endTime = str(dt.datetime.now())
            
class FunctionStep(Step):
    """This is a Step wrapper around a normal function"""
    def __init__(self, func=None, *funcArgs):
        """Receive the function to execute and the 
        parameters to call it"""
        Step.__init__(self)
        self.func = func
        self.funcArgs = funcArgs
        self.funcName = String(func)
        self.argsStr = String(pickle.dumps(funcArgs))
        
    def _run(self):
        self.func(*self.funcArgs)
        

#def ProtocolType(type):
#    """Protocols metaclass"""
#    def __init__(cls, name, bases, dct):
#        print '-----------------------------------'
#        print "Initializing protocol class", name
#        print dct
#        super(ProtocolType, cls).__init__(name, bases, dct)
        
def loadDef():
    print "loading..."
                
class Protocol(Step):
    """The Protocol is a higher type of Step.
    It also have the inputs, outputs and other Steps properties,
    but contains a list of steps that are executed"""
    #__metaclass__ = ProtocolType
    # Params definition for this class
    _paramDefinition = loadDef()
    
    def __init__(self, **args):
        Step.__init__(self, **args)
        self.addAttribute('steps', List)
        self.workingDir = args.get('workingDir', '.')
        self.mapper = args.get('mapper', None)
        
    def store(self, obj, commit=False):
        if not self.mapper is None:
            self.mapper.store(obj)
            if commit:
                self.commit()
            
    def commit(self):
        if not self.mapper is None:
            self.mapper.commit()
            
    def defineSteps(self):
        """Define all the steps that will be executed."""
        pass
    
    def insertStep(self, step):
        """Insert a new step in the list"""
        self.steps.append(step)
        
    def insertFunctionStep(self, func, *funcArgs):
        self.insertStep(FunctionStep(func, *funcArgs))
        
    def _run(self):
        self.defineSteps()
        
        for step in self.steps:
            step.run()
            print "step_index: ", self.steps.getIndexStr(self.currentStep)
            self.mapper.insertChild(self.steps, self.steps.getIndexStr(self.currentStep),
                             step, self.namePrefix)
            self.currentStep += 1
            self.commit()
            if step.status == STATUS_FAILED:
                raise Exception("Protocol failed")

    def run(self):
        self.currentStep = 1
        self.store(self, commit=True)
        print "self.steps.name: ", self.steps.name
        self.namePrefix = replaceExt(self.steps.name, str(self.steps.id)) #keep 
        Step.run(self)
        print 'PROTOCOL FINISHED'
        #self.store(self, commit=True)

class ProtImportMicrographs(Protocol):
    def __init__(self, **args):
        Protocol.__init__(self, **args)
        self.addAttribute('pattern', String)
        self.addAttribute('n', Integer)
        self.count = 1
        
    def defineSteps(self):
        for i in range(self.n):
            self.insertFunctionStep(self.copyMicrographs, self.pattern)
        
    def copyMicrographs(self, pattern):
        """Copy micrograph with filename matching the pattern"""
        import time
        print "step: ", self.count
        time.sleep(1)
        self.count += 1
        

class ProtScreenMicrographs(Protocol):
    pass

class ProtDownsampleMicrographs(Protocol):
    pass

class ProtParticlePicking(Protocol):
    pass

class ProtAlign(Protocol):
    pass

class ProtClassify(Protocol):
    pass

class ProtAlignClassify(Protocol):
    pass
