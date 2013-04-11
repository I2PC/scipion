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


LEVEL_NORMAL = 0
LEVEL_ADVANCED = 1
LEVEL_EXPERT = 2


class FormBase(OrderedObject):
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
        self.label = String(args.get('label', None))
        self.expertLevel = Integer(args.get('expertLevel', LEVEL_NORMAL))
        
    def isExpert(self):
        return self.expert.hasValue()
    
    
class Section(FormBase):
    """Definition of a section to hold other params"""
    def __init__(self, **args):
        FormBase.__init__(self, **args)
        self.questionParam = String(args.get('questionParam', None))
        self.paramList = []
    
    def addParam(self, name, ParamClass, **args):
        """Add a new param to last section"""
        p = ParamClass(**args)
        setattr(self, name, p)
        self.paramList.append(name)
        return p
    
    def getChildParams(self):
        """Return key and param for every child param"""
        for name in self.paramList:
            yield (name, getattr(self, name))
        
    def __str__(self):
        s = "  Section, label: %s\n" % self.label.get()
        for key in self._attributes:
            s += "    Param: '%s', %s\n" % (key, str(getattr(self, key)))
        return s        

                    
class Form(List):
    """Store all sections and parameters"""

    def addSection(self, **args):
        """Add a new section"""
        s = Section(**args)
        self.append(s)
        return s

    def addParam(self, name, ParamClass, **args):
        """Add a new param to last section"""
        return self[-1].addParam(name, ParamClass, **args)
        
    def __str__(self):
        s = "Form: \n"
        for section in self:
            s += str(section)
        return s
        
      
class Param(FormBase):
    """Definition of a protocol parameter"""
    def __init__(self, **args):
        FormBase.__init__(self, **args)
        self.paramClass = args.get('paramClass', None) # This should be defined in subclasses
        self.default = String(args.get('default', None))
        self.help = String(args.get('help', None))
        self.isImportant = args.get('important', False)
        
    def __str__(self):
        return "    label: %s" % self.label.get()


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
        self.pointerClass = args.get('pointerClass')
        
        
        
        
        
