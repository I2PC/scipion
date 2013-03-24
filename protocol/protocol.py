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

from pyworkflow.object import FakedObject, String, List

STATUS_LAUNCHED = "launched"  # launched to queue system
STATUS_RUNNING = "running"    # currently executing
STATUS_FAILED = "failed"      # it have been failed
STATUS_FINISHED = "finished"  # successfully finished
STATUS_WAITING = "waiting"    # waiting for user interaction

class Step(FakedObject):
    """Basic execution unit.
    It should defines its Input, Output
    and define a run method"""
    def __init__(self):
        FakedObject.__init__(self)
        self._inputs = []
        self._outputs = []
        self.addAttribute('status', String)
        self.addAttribute('initTime', String)
        self.addAttribute('endTime', String)
        
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
        except Exception, e:
            self.status = STATUS_FAILED
        finally:
            self.endTime = str(dt.datetime.now())
            

def ProtocolType(type):
    """Protocols metaclass"""
    def __init__(cls, name, bases, dct):
        print '-----------------------------------'
        print "Initializing protocol class", name
        print dct
        super(ProtocolType, cls).__init__(name, bases, dct)
        
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
        self._steps = List() # The list of steps
        self.workingDir = args.get('workingDir', None)
        
    def insertStep(self, step):
        """Insert a new step in the list"""
        self._steps.append(step)
        
    def defineSteps(self):
        """Define all the steps that will be executed."""
        pass
    
    def run(self):
        for step in self._steps:
            step.run()
    

class ProtImportMicrographs(Protocol):
    pass

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
