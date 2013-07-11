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
This module have the classes for protocol params definition:
Param, Section and Form
The definition will be holded at class level and will
be shared by all protocol class instances
"""

from pyworkflow.object import *
import re
import collections


LEVEL_NORMAL = 0
LEVEL_ADVANCED = 1
LEVEL_EXPERT = 2
LEVEL_CHOICES = ('Normal', 'Advanced', 'Expert')


class FormElement(OrderedObject):
    """Base for any element on the form"""
    
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
        self.label = String(args.get('label', None))
        self.expertLevel = Integer(args.get('expertLevel', LEVEL_NORMAL))
        self.condition = String(args.get('condition', None))
        # This two list will be filled by the Form
        # which have a global-view of all parameters
        self._dependants = [] # All param names in which condition appears this param
        self._conditionParams = [] # All param names that appears in the condition
        
    def isExpert(self):
        return self.expert.hasValue()
    
    def hasCondition(self):
        return self.condition.hasValue()
    
        
class Param(FormElement):
    """Definition of a protocol parameter"""
    def __init__(self, **args):
        FormElement.__init__(self, **args)
        self.paramClass = args.get('paramClass', None) # This should be defined in subclasses
        self.default = String(args.get('default', None))
        self.help = String(args.get('help', None))
        self.isImportant = Boolean(args.get('important', False))
        self.validators = args.get('validators', [])
        
    def __str__(self):
        return "    label: %s" % self.label.get()
    
    def addValidator(self, validator):
        """ Validators should be callables that 
        receive a value and return a list of errors if so.
        If everything is ok, the result should be an empty list.
        """
        self.validators.append(validator)
        
    def validate(self, value):
        errors = []
        for val in self.validators:
            errors += val(value)
        return errors
    
    
class Section(FormElement):
    """Definition of a section to hold other params"""
    def __init__(self, form, **args):
        FormElement.__init__(self, **args)
        self._form = form
        self.questionParam = String(args.get('questionParam', ''))
        self._paramList = []
    
    def hasQuestion(self):
        """Return True if a question param was set"""
        return self.questionParam.get() in self._paramList
    
    def getQuestionName(self):
        """ Return the name of the question param. """
        return self.questionParam.get()
    
    def getQuestion(self):
        """Return the question param"""
        return self._form.getParam(self.questionParam.get())
    
    def addParam(self, name):
        """Add a new param to last section"""
        self._paramList.append(name)
    
    def iterParams(self):
        """Return key and param for every child param"""
        for name in self._paramList:
            yield (name, self._form.getParam(name))

                    
class Form():
    """Store all sections and parameters"""
    def __init__(self):
        self._sectionList = [] # Store list of sections
        self._paramsDict = collections.OrderedDict() #{} # Dictionary to store all params, grouped by sections
        self._lastSection = None
        self.addGeneralSection()
        
    def addSection(self, **args):
        """Add a new section"""
        self.lastSection = Section(self, **args)
        self._sectionList.append(self.lastSection)
        return self.lastSection

    def addParam(self, paramName, ParamClass, **args):
        """Add a new param to last section"""
        param = ParamClass(**args)
        self._paramsDict[paramName] = param
        self._analizeCondition(paramName, param)
        return self.lastSection.addParam(paramName)
    
    def _analizeCondition(self, paramName, param):
        if param.hasCondition():
            param._conditionParams = []
            tokens = re.split('\W+', param.condition.get())
            for t in tokens:
                if self.hasParam(t):
                    param._conditionParams.append(t)
                    self.getParam(t)._dependants.append(paramName)
                    #print "tokens found: ", self.getParam(t)._dependants
#                else:
#                    raise Exception("Invalid token '%s' used in param '%s' condition" 
#                                    % (t, paramName))
    def evalCondition(self, protocol, paramName):
        """Evaluate if a condition is True for a give param
        with the values of a particular Protocol"""
        param = self.getParam(paramName)
        if not param.hasCondition():
            return True
        condStr = param.condition.get()
        for t in param._conditionParams:
            if self.hasParam(t):
                condStr = condStr.replace(t, str(protocol.getAttributeValue(t)))
        return eval(condStr)
    
    def validateParams(self, protocol):
        """ Check that all validations of the params in the form
        are met for the protocol param values.
        It will return a list with errors, just in the same
        way of the Protocol.validate function
        """
        errors = []
        
        for name, param in self.iterParams():
            value = protocol.getAttributeValue(name)
            errors += param.validate(value)
        
        return errors
        
        
    def getParam(self, paramName):
        """Retrieve a param given a the param name
        None is returned if not found
        """        
        return self._paramsDict.get(paramName, None)
    
    def hasParam(self, paramName):
        return paramName in self._paramsDict
        
    def __str__(self):
        s = "Form: \n"
        for section in self.iterSections():
            s += str(section)
        return s
    
    def iterSections(self):
        return self._sectionList
    
    def iterParams(self):
        return self._paramsDict.iteritems()
    
    def iterPointerParams(self):
        for paramName, param in self._paramsDict.iteritems():
            if isinstance(param, PointerParam):
                yield paramName, param

    def addGeneralSection(self):
        self.addSection(label='General')
        self.addParam('runName', StringParam, label="Run name:", important=True, 
                      help='Select run name label to identify this run.')
        self.addParam('showComment', BooleanParam, default=False, 
                      label="Show comment?")
        self.addParam('comment', StringParam, condition="showComment",
                      label="Comment:", help='Make some annotations on this run.')
        self.addParam('runMode', EnumParam, choices=['resume', 'restart'],
                      label="Run mode", display=EnumParam.DISPLAY_COMBO, default=0,
                      help='The <resume> mode will try to start the execution'
                           'from the last sucessfully finished step if possible.'
                           'On the contrary, <restart> will delete all previous results'
                           'of this particular run and start from the beginning. This option'
                           'should be used carefully.'
                      )
  
    def addParallelSection(self, threads=1, mpi=8, condition="",
                           hours=72, jobsize=0):
        self.addSection(label='Parallelization')
        self.addParam('hostName', StringParam, default="localhost",
                      label='Execution host',
                      help='Select in which of the available do you want to launch this protocol.')
        if threads > 0:
            self.addParam('numberOfThreads', IntParam, default=threads,
                          label='Threads',
                          help='This option provides shared-memory parallelization on multi-core machines.'
                                'It does not require any additional software, other than <Xmipp>' )
        if mpi > 0:
            self.addParam('numberOfMpi', IntParam, default=mpi,
                          label='MPI processes',
                          help='This option provides the number of independent processes spawned'
                                'in parallel by <mpirun> command in a cluster, usually throught'
                                'a queue system. This will require that you have compile <Xmipp>'
                                'with <mpi> support.')
        if jobsize > 0:
            self.addParam('mpiJobSize', IntParam, default=jobsize,
                          label='MPI job size', condition="numberOfMpi>1",
                          help='Minimum size of jobs in mpi processes.'
                               'Set to 1 for large images (e.g. 500x500)'
                               'and to 10 for small images (e.g. 100x100)')

        

# More Param sub-classes

class StringParam(Param):
    """Param with underlying String value"""
    def __init__(self, **args):
        Param.__init__(self, paramClass=String, **args)


class TextParam(StringParam):
    """Long string params"""
    def __init__(self, **args):
        StringParam.__init__(self, **args)
        
        
class RegexParam(StringParam):
    """Regex based string param"""
    pass


class PathParam(StringParam):
    """Param for path strings"""
    pass


class FileParam(PathParam):
    """Filename path"""
    pass


class FolderParam(PathParam):
    """Folder path"""
    pass

        
class IntParam(Param):
    def __init__(self, **args):
        Param.__init__(self, paramClass=Integer, **args)
        
        
class EnumParam(IntParam):
    """Select from a list of values, separated by comma"""
    # Possible values for display
    DISPLAY_LIST = 0
    DISPLAY_COMBO = 1
    
    def __init__(self, **args):
        IntParam.__init__(self, **args)
        self.choices = args.get('choices', [])
        self.display = Integer(args.get('display', EnumParam.DISPLAY_LIST))
    
    
class FloatParam(Param):
    def __init__(self, **args):
        Param.__init__(self, paramClass=Float, **args)
        
        
class BooleanParam(Param):
    def __init__(self, **args):
        Param.__init__(self, paramClass=Boolean, **args)

        
class PointerParam(Param):
    def __init__(self, **args):
        Param.__init__(self, paramClass=Pointer, **args)
        self.pointerClass = String(args.get('pointerClass'))
        
        

# ------------------------------------------------------------------------
#         Validators
#-------------------------------------------------------------------------

class Validator(object):
    pass

class Conditional(Validator):
    """ Simple validation based on a condition. 
    If the value doesn't meet the condition,
    the error will be returned.
    """
    def __init__(self, error):
        self.error = error
        
    def __call__(self, value):
        errors = []
        if not self._condition(value):
            errors.append(self.error)
        return errors    
    
class LT(Conditional):
    def __init__(self, thresold, error='Value should be less than the thresold'):
        Conditional.__init__(self, error)
        self._condition = lambda value: value < thresold
        
class LE(Conditional):
    def __init__(self, thresold, error='Value should be less or equal than the thresold'):
        Conditional.__init__(self, error)
        self._condition = lambda value: value <= thresold        
        
class GT(Conditional):
    def __init__(self, thresold, error='Value should be greater than the thresold'):
        Conditional.__init__(self, error)
        self._condition = lambda value: value > thresold

class GE(Conditional):
    def __init__(self, thresold, error='Value should be greater or equal than the thresold'):
        Conditional.__init__(self, error)
        self._condition = lambda value: value >= thresold               

class Range(Conditional):
    def __init__(self, minValue, maxValue, error='Value is outside range'):
        Conditional.__init__(self, error)
        self._condition = lambda value: value >= minValue and value <= maxValue
        
class Positive(GT):
    def __init__(self, error='Value should be greater than zero'):
        GT.__init__(self, 0.0, error)

            

        
        
