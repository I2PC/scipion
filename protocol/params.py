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


class FormBase(OrderedObject):
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
         
    def setTags(self, tags):
        for k, v in tags.iteritems():
            if k in ['section']:
                value = None
            elif k in ['hidden', 'expert', 'view', 'has_question', 'visualize', 'text', 'file']:
                value = Boolean(True)
            else:
                value = String(v)
            if not value is None:
                setattr(self, k, value)
                
    def __getattr__(self, name):
        value = None
        if name in ['hidden', 'expert', 'view', 'has_question', 'visualize', 'text', 'file']:
            value = Boolean(True)
        elif name in ['label', 'help', 'list', 'condition', 'validate', 'wizard']:
            value = String()
        elif name == 'default':
            value = DefaultString()
            
        self.__setattr__(name, value)
        return value
    
      
class Param(FormBase):
    """Definition of a protocol paramter"""
    pass
        
    
class Section(FormBase):
    """Definition of a section to hold other params"""
    
    def addParam(self, name, value):
        """Add a new param to last section"""
        setattr(self, name, value)        

                    
class Form(List):
    """Store all sections and parameters"""

    def addSection(self, section):
        """Add a new section"""
        self.lastSection = section
        List.append(self, section)

    def addParam(self, name, value):
        """Add a new param to last section"""
        self.lastSection.addParam(name, value)

class DefaultString(String):
    pass 
        
        
        
        
        
        