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
This module contains reflection utilities
(dynamically load classes, inspect object properties and others) 
"""

import os, sys
from os.path import exists, join
from inspect import isclass
    
def getModules(path):
    """ Try to find possible sub-modules under path.
    A dictionary will be returned with modules names
    as keys and the modules objects as values.
    """
    sys.path.append(path)
    folders = os.listdir(path)
    modules = {}
    
    for f in folders:
        if exists(join(path, f, '__init__.py')):
            try:
                m = __import__(f)
                modules[f] = m
            except Exception:
                print "Error loading module: '%s'" % f
    
    return modules
    
def getSubclassesFromModules(BaseClass, modules):
    """ Find subclasses of BaseClass from a give dict of modules.
    """
    subclasses = {}
    
    for m in modules.values():
        subDict = getSubclasses(BaseClass, m.__dict__)
        subclasses.update(subDict)
    
    return subclasses

def getSubclassesFromPath(BaseClass, path):
    """ Try to find possible sub-packages under path
    and find subclasses of BaseClass from them
    Return a dictionary containing the subclasses.
    """
    modules = getModules(path)
    return getSubclassesFromModules(BaseClass, modules)
    
def getSubclasses(BaseClass, inputDict):
    """ Iterate over inputDict and find all subclasses
    of BaseClass, that will be set in outputDict.
    """
    outputDict = {}
    for k, v in inputDict.iteritems():
        if isclass(v) and issubclass(v, BaseClass):
            outputDict[k] = v
    return outputDict


    