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

class Param(Object):
    """Definition of a protocol paramter"""
    
    def __init__(self, **args):
        """Create a new Param. Possible values in **args:
        paramClass: this will be the class to store the param value(Integer, String...)
        paramLabel: this will be used to display a text about the param
           helpMsg: string value with the help about the param
              tags: a dictionary containing several tags
        conditions: some dependencies with other params
        validators: some validations imposed to the param"""
        Object.__init__(self, **args)
    
class Section(Object):
    """Definition of a section to hold other params"""
    
    def __init__(self, **args):
        """Posible values in **args are:
        sectionLabel: string for the section label
             helpMsg: help string message"""
        Object.__init__(self, **args)
             
        
        
        
        
        
        