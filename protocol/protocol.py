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

import datetime as dt
import pickle

from pyworkflow.object import Object, String, List, Integer
from pyworkflow.utils.path import replaceExt

STATUS_LAUNCHED = "launched"  # launched to queue system
STATUS_RUNNING = "running"    # currently executing
STATUS_FAILED = "failed"      # it have been failed
STATUS_FINISHED = "finished"  # successfully finished
STATUS_WAITING = "waiting"    # waiting for user interaction

class Step(Object):
    """Basic execution unit.
    It should defines its Input, Output
    and define a run method"""
    def __init__(self, **args):
        Object.__init__(self, **args)
        self._inputs = []
        self._outputs = []
        self.status = String()
        self.initTime = String()
        self.endTime = String()
        self.error = String()
        
    def _storeAttributes(self, attrList, attrDict):
        """Store all attributes in attrDict as 
        attributes of self, also store the key in attrList"""
        for key, value in attrDict.iteritems():
            attrList.append(key)
            setattr(self, key, value)
            
    def _getAttribute(self, name):
        return self._attributes.get(name, None)
        
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
        self.initTime.set(dt.datetime.now())
        self.endTime.set(None)
        try:
            self._run()
            self.status.set(STATUS_FINISHED)
        except Exception, e:
            self.status.set(STATUS_FAILED)
            self.error.set(e)
            raise #only in development
        finally:
            self.endTime.set(dt.datetime.now())
            
class FunctionStep(Step):
    """This is a Step wrapper around a normal function
    This class will ease the insertion of Protocol function steps
    throught the function insertFunctionStep"""
    def __init__(self, func=None, *funcArgs):
        """Receive the function to execute and the 
        parameters to call it"""
        Step.__init__(self)
        self.func = func
        self.funcArgs = funcArgs
        self.funcName = String(func)
        self.argsStr = String(pickle.dumps(funcArgs))
        print "funcName", self.funcName.get()
        print "funcArgs:", funcArgs
        print "args: ", pickle.loads(self.argsStr.get())
        
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
        self.steps = List()
        self.workingDir = args.get('workingDir', '.')
        self.mapper = args.get('mapper', None)
        
    def store(self, *objs):
        if not self.mapper is None:
            if len(objs) == 0:
                self.mapper.store(self)
            else:
                for obj in objs:
                    print "storing: ", obj
                    self.mapper.store(obj)
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
            self.mapper.insertChild(self.steps, self.steps.getIndexStr(self.currentStep),
                             step, self.namePrefix)
            self.currentStep += 1
            self.mapper.commit()
            if step.status == STATUS_FAILED:
                raise Exception("Protocol failed")

    def run(self):
        self.currentStep = 1
        self.store()
        self.namePrefix = replaceExt(self.steps.getName(), self.steps.strId()) #keep 
        Step.run(self)
        self.store(self.status, self.initTime, self.endTime)
        print 'PROTOCOL FINISHED'

class ProtImportMicrographs(Protocol):
    def __init__(self, **args):
        Protocol.__init__(self, **args)
        self.pattern = String()
        self.n = Integer()
        self.count = 1
        
    def defineSteps(self):
        for _ in range(self.n):
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
    
    def __init__(self, **args):
        Protocol.__init__(self, **args)
    

class ProtAlign(Protocol):
    pass

class ProtClassify(Protocol):
    pass

class ProtAlignClassify(Protocol):
    pass
