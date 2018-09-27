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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains reflection utilities
(dynamically load classes, inspect object properties and others)
"""

import os
from os.path import exists, join
import sys
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
                checkPlugin(m)
            except Exception as ex:
                print(">>> Error loading module: '%s'" % f)
                print(">>> Exception: ", ex)
                import traceback
                traceback.print_exc()

    return modules


def getSubclassesFromModules(BaseClass, modules, debug=False):
    """ Find subclasses of BaseClass from a give dict of modules.
    """
    subclasses = {}

    for m in modules.values():
        if debug:
            print("loading module: ", m.__name__)
        subDict = getSubclasses(BaseClass, m.__dict__)
        
        for subclass in subDict.values():
            # some protocols have pyworkflow.em.packages. in the __module__ and other no
            moduleName = subclass.__module__.replace('pyworkflow.em.packages.', '')
            if moduleName.startswith(m.__name__):
                subclass._package = m
                if debug:
                    print("  found: ", subclass.__name__, "module: ", subclass.__module__)
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


def checkPlugin(module):
    if not getattr(module, '_plugin', None):
        print('WARNING: module "%s" using old package structure, '
              '_plugin attribute missing' % module.__name__)