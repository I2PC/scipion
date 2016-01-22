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

import re
import collections

from pyworkflow.object import *
from constants import *


class FormElement(OrderedObject):
    """Base for any element on the form"""
    ATTRIBUTES = ['label', 'expertLevel', 'condition', 'important', 'help',
                  'default', 'paramClass']
    
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
        self.label = String(args.get('label', None))
        self.expertLevel = Integer(args.get('expertLevel', LEVEL_NORMAL))
        self.condition = String(args.get('condition', None))
        self._isImportant = Boolean(args.get('important', False))
        self.help = String(args.get('help', None))
        # This two list will be filled by the Form
        # which have a global-view of all parameters
        self._dependants = [] # All param names in which condition appears this param
        self._conditionParams = [] # All param names that appears in the condition
        
    def isExpert(self):
        return self.expertLevel > LEVEL_NORMAL
    
    def setExpert(self):
        self.expertLevel.set(LEVEL_ADVANCED)
        
    def isImportant(self):
        return self._isImportant.get()
    
    def setImportant(self, value):
        self._isImportant.set(value)
    
    def hasCondition(self):
        return self.condition.hasValue()
    
    def getLabel(self):
        return self.label.get()
    
    def config(self, **kwargs):
        """ Configure the object and set attributes
        coming in the keyword-arguments, the 
        same as in the __init__
        """
        for key in self.ATTRIBUTES:
            if key in kwargs:
                self.setAttributeValue(key, kwargs.get(key)) 
    
        
class Param(FormElement):
    """Definition of a protocol parameter"""
    def __init__(self, **args):
        FormElement.__init__(self, **args)
        self.paramClass = args.get('paramClass', None) # This should be defined in subclasses
        self.default = String(args.get('default', None))
        #from pyworkflow.web.app.views_util import parseText
        #self.help = parseText(String(args.get('help', None)))

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
    
    def getDefault(self):
        return self.default.get()
    
    def setDefault(self, newDefault):
        self.default.set(newDefault)
    
    
class ElementGroup(FormElement):
    """ Class to group some params in the form.
    Such as: Labeled group or params in the same line. 
    """
    def __init__(self, form=None, **args):
        FormElement.__init__(self, **args)
        self._form = form
        self._paramList = []    
    
    def iterParams(self):
        """ Return key and param for every child param. """
        for name in self._paramList:
            yield (name, self._form.getParam(name))

    def addParam(self, paramName, ParamClass, **kwargs):
        """Add a new param to the group"""
        param = ParamClass(**kwargs)
        self._paramList.append(paramName)
        self._form.registerParam(paramName, param)
        return param
    
    def addHidden(self, paramName, ParamClass, **kwargs):
        """Add a hidden parameter to be used in conditions. """
        kwargs.update({'label': '', 'condition': 'False'})
        self.addParam(paramName, ParamClass, **kwargs)

    def addLine(self, lineName, **kwargs):
        
        labelName = lineName
        for symbol in ' ()':
            labelName = labelName.replace(symbol, '_')
        
        return self.addParam(labelName, Line, form=self._form, 
                             label=lineName, **kwargs)        
    
    
# ----------- Some type of ElementGroup --------------------------

class Line(ElementGroup):
    """ Group to put some parameters in the same line. """
    pass


class Group(ElementGroup):
    """ Group some parameters with a labeled frame. """
    pass

    
class Section(ElementGroup):
    """Definition of a section to hold other params"""
    def __init__(self, form, **args):
        ElementGroup.__init__(self, form, **args)
        self.questionParam = String(args.get('questionParam', ''))
    
    def hasQuestion(self):
        """Return True if a question param was set"""
        return self.questionParam.get() in self._paramList
    
    def getQuestionName(self):
        """ Return the name of the question param. """
        return self.questionParam.get()
    
    def getQuestion(self):
        """ Return the question param"""
        return self._form.getParam(self.questionParam.get())

    def addGroup(self, groupName, **kwargs):
        labelName = groupName
        for symbol in ' ()':
            labelName = labelName.replace(symbol, '_')
        
        return self.addParam(labelName, Group, form=self._form, 
                             label=groupName, **kwargs)
        
class Form(object):
    """Store all sections and parameters"""
    def __init__(self, protocol):
        """ Build a Form from a given protocol. """
        object.__init__(self)
        self._sectionList = [] # Store list of sections
        self._paramsDict = collections.OrderedDict() #{} # Dictionary to store all params, grouped by sections
        self._lastSection = None
        self._protocol = protocol
        self.addGeneralSection()
        
    def getClass(self):
        return type(self)
        
    def addSection(self, label='', **kwargs):
        """Add a new section"""
        self.lastSection = Section(self, label=label, **kwargs)
        self._sectionList.append(self.lastSection)
        return self.lastSection
    
    def addGroup(self, *args, **kwargs):
        return self.lastSection.addGroup(*args, **kwargs)
    
    def addLine(self, *args, **kwargs):
        return self.lastSection.addLine(*args, **kwargs)      

    def registerParam(self, paramName, param):
        """ Register a given param in the form. """
        self._paramsDict[paramName] = param        
        self._analizeCondition(paramName, param)
        
    def addParam(self, *args, **kwargs):
        """Add a new param to last section"""
        return self.lastSection.addParam(*args, **kwargs)
    
    def addHidden(self, *args, **kwargs):
        return self.lastSection.addHidden(*args, **kwargs)
    
    def _analizeCondition(self, paramName, param):
        if param.hasCondition():
            param._conditionParams = []
            tokens = re.split('\W+', param.condition.get())
            for t in tokens:
                if self.hasParam(t):
                    self.getParam(t)._dependants.append(paramName)
                    param._conditionParams.append(t)
                if self._protocol.hasAttribute(t):
                    param._conditionParams.append(t)
                    
    def escapeLiteral(self, value):
        if isinstance(value, str):
            result = "'%s'" % value
        else:
            result = str(value)
        return result
    
    def evalParamCondition(self, paramName):
        """Evaluate if a condition is True for a give param
        with the values of a particular Protocol"""
        param = self.getParam(paramName)
        if not param.hasCondition():
            return True
        condStr = param.condition.get()
        for t in param._conditionParams:
            if self.hasParam(t) or self._protocol.hasAttribute(t):
                condStr = condStr.replace(t, self.escapeLiteral(self._protocol.getAttributeValue(t)))
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
    
    def iterAllParams(self):
        """ Iter all parameters, including ElementGroups. """
        return self._paramsDict.iteritems()
    
    def iterParams(self):
        """ Iter parameters disregarding the ElementGroups. """
        for k, v in self._paramsDict.iteritems():
            if not isinstance(v, ElementGroup):
                yield k, v
        
    def iterPointerParams(self):
        for paramName, param in self._paramsDict.iteritems():
            if isinstance(param, PointerParam):
                yield paramName, param

    def addGeneralSection(self):
        self.addSection(label='General')
        self.addParam('runName', StringParam, label="Run name:", important=True, 
                      help='Select run name label to identify this run.')
        self.addParam('runMode', EnumParam, choices=['resume', 'restart'],
                      label="Run mode", display=EnumParam.DISPLAY_COMBO, default=0,
                      help='The <resume> mode will try to start the execution'
                           'from the last successfully finished step if possible.'
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

# TODO: Handle filter pattern
class FileParam(PathParam):
    """Filename path"""
    pass


class FolderParam(PathParam):
    """Folder path"""
    pass


class LabelParam(StringParam):
    """ Just the same as StringParam, but to be rendered
    as a label and can not be directly edited by the user
    in the Protocol Form.
    """
    pass

        
class IntParam(Param):
    def __init__(self, **args):
        Param.__init__(self, paramClass=Integer, **args)
        self.addValidator(Format(int, error="should be an integer",
                                 allowsNull=args.get('allowsNull', False)))
        
        
class EnumParam(IntParam):
    """Select from a list of values, separated by comma"""
    # Possible values for display
    DISPLAY_LIST = 0
    DISPLAY_COMBO = 1
    DISPLAY_HLIST = 2 # horizontal list, save space
    
    def __init__(self, **args):
        IntParam.__init__(self, **args)
        self.choices = args.get('choices', [])
        self.display = Integer(args.get('display', EnumParam.DISPLAY_COMBO))
    
    
class FloatParam(Param):
    def __init__(self, **args):
        Param.__init__(self, paramClass=Float, **args)
        self.addValidator(Format(float, error="should be a float",
                                 allowsNull=args.get('allowsNull', False)))

        
class BooleanParam(Param):
    def __init__(self, **args):
        Param.__init__(self, paramClass=Boolean, **args)
        self.addValidator(NonEmptyBool)


class HiddenBooleanParam(BooleanParam):
    def __init__(self, **args):
        Param.__init__(self, paramClass=Boolean, **args)

        
class PointerParam(Param):
    """ This type of Param will serve to select existing objects
    in the database that will be input for some protocol.
    """
    def __init__(self,  paramClass=Pointer, **args):
        Param.__init__(self, paramClass=paramClass, **args)
        # This will be the class to be pointed
        pointerClass = args.get('pointerClass')
        if ',' in pointerClass:
            self.pointerClass = CsvList()
            self.pointerClass.set(pointerClass)
        else:
            self.pointerClass = String(pointerClass) 

        # Some conditions on the pointed candidates
        self.pointerCondition = String(args.get('pointerCondition', None))
        self.allowsNull = Boolean(args.get('allowsNull', False))
        
    def setPointerClass(self, newPointerClass):
        self.pointerClass.set(newPointerClass)


class MultiPointerParam(PointerParam):
    """ This type of Param will serve to select objects
    with DIFFERENT types from the database to be input for some protocol.
    """
    def __init__(self, **args):
        PointerParam.__init__(self, paramClass=PointerList, **args)
        self.maxNumObjects = Integer(args.get('maxNumObjects',100))
        self.minNumObjects = Integer(args.get('minNumObjects',2))   

        
class RelationParam(Param):
    """ This type of Param is very similar to PointerParam, since it will
    hold a pointer to another object. But, in the PointerParam, we search
    for objects of some Class (maybe with some conditions).
    Here, we search for objects related to a given attribute of a protocol
    by a given relation.
    """
    def __init__(self, **args):
        Param.__init__(self, paramClass=Pointer, **args)
        # This will be the name of the relation
        self._relationName = String(args.get('relationName'))
        # We will store the attribute name in the protocol to be 
        # used as the object for which relations will be search
        self._attributeName = String(args.get('attributeName'))
        # This specify if we want to search for childs or parents
        # of the given attribute of the protocol
        self._direction = Integer(args.get('direction', RELATION_CHILDS))
        self.allowsNull = Boolean(args.get('allowsNull', False))
        
    def getName(self):
        return self._relationName.get()
    
    def getAttributeName(self):
        return self._attributeName.get()
    
    def getDirection(self):
        return self._direction.get()       
        
        
class ProtocolClassParam(StringParam):
    def __init__(self, **args):
        StringParam.__init__(self, **args)
        self.protocolClassName = String(args.get('protocolClassName'))
        self.allowSubclasses = Boolean(args.get('allowSubclasses', False))
        
        
class DigFreqParam(FloatParam):
    """ Digital frequency param. """
    def __init__(self, **args):
        FloatParam.__init__(self, **args)
        self.addValidator(FreqValidator)
        
        
class NumericListParam(StringParam):
    """ This class will serve to have list representations as strings.
     Possible notation are:
     1000 10 1 1 -> to define a list with 4 values [1000, 10, 1, 1], or
     10x2 5x3    -> to define a list with 5 values [10, 10, 5, 5, 5]
     If you ask for more elements than in the list, the last one is repeated
    """
    def __init__(self, **args):
        StringParam.__init__(self, **args)
        self.addValidator(NumericListValidator())
        
        
class NumericRangeParam(StringParam):
    """ This class will serve to specify range of numbers with a string representation.
     Possible notation are:
        "1,5-8,10" -> [1,5,6,7,8,10]
        "2,6,9-11" -> [2,6,9,10,11]
        "2 5, 6-8" -> [2,5,6,7,8]
    """
    def __init__(self, **args):
        StringParam.__init__(self, **args)
        # TODO: ADD a syntax validator
        
        
class TupleParam(Param):
    """ This class will condense a tuple of several
    other params of the same type and related concenpt.
    For example: min and max, low and high.
    """
    def __init__(self, **args):
        Param.__init__(self, **args)


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
    def __init__(self, error, allowsNull=False):
        self.error = error
        self._allowsNull = allowsNull
        
    def __call__(self, value):
        errors = []
        if (value is not None or not self._allowsNull):
            if not self._condition(value):
                errors.append(self.error)
        return errors   
    
    
class Format(Conditional):
    """ Check if the format is right. """
    def __init__(self, valueType, error='Value have not a correct format', allowsNull=False):
        Conditional.__init__(self, error, allowsNull)
        self.valueType = valueType
        
    def _condition(self, value):
        try:
            self.valueType(value)
            return True
        except Exception:
            return False


class NonEmptyCondition(Conditional):
    def __init__(self, error='Value cannot be empty'):
        Conditional.__init__(self, error)
        self._condition = lambda value: len(value) > 0
        
        
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
        
        
class NumericListValidator(Conditional):
    """ Validator for ListParam. See ListParam. """
    def __init__(self, error='Incorrect format for numeric list param. '):
        Conditional.__init__(self, error)
        
    def _condition(self, value):
        try:
            value = value.replace('x', '')
            parts = value.split()
            for p in parts:
                float(p)
            return True
        except Exception:
            return False    


class NonEmptyBoolCondition(Conditional):
    def __init__(self, error='Boolean param needs to be set.'):
        Conditional.__init__(self, error)
        self._condition = lambda value: value is not None

#--------- Some constants validators ---------------------

Positive = GT(0.0, error='Value should be greater than zero')

FreqValidator = Range(0., 0.5, 
                      error="Digital frequencies should be between 0. and 0.5")

NonEmpty = NonEmptyCondition()
NonEmptyBool = NonEmptyBoolCondition()
