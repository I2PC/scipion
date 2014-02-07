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
This modules contains classes related with EM
"""
import os

import pyworkflow as pw
from pyworkflow.utils.reflection import getSubclassesFromModules, getSubclasses, getModules
from data import *
from protocol import *
from constants import *
from convert import *
from pyworkflow.viewer import Viewer
from pyworkflow.wizard import Wizard
from viewer import *
#from packages import *

PACKAGES_PATH = os.path.join(pw.HOME, 'em', 'packages')
PACKAGES_DICT = getModules(PACKAGES_PATH)

# Load all Protocol subclasses found in EM-packages
emProtocolsDict = getSubclassesFromModules(Protocol, PACKAGES_DICT)
emProtocolsDict.update(getSubclasses(Protocol, globals()))

# Load all EMObject subclasses found in EM-packages
emObjectsDict = getSubclassesFromModules(EMObject, PACKAGES_DICT)
emObjectsDict.update(getSubclasses(EMObject, globals()))

# Load all subclasses of Viewer of different packages
emViewersDict = getSubclassesFromModules(Viewer, PACKAGES_DICT)

# Load all subclasses of Wizards
emWizardsDict = getSubclassesFromModules(Wizard, PACKAGES_DICT)
        
def findClass(className):
    if className in emProtocolsDict:
        return emProtocolsDict[className]
    if className in emObjectsDict:
        return emObjectsDict[className]
    raise Exception("findClass: class '%s' not found." % className)

def findSubClasses(classDict, className):
    """ Find all subclasses of a give className. """
    cls = classDict[className]
    subclasses = {}
    
    for k, v in classDict.iteritems():
        if issubclass(v, cls):
            subclasses[k] = v    
    return subclasses


def findViewers(className, environment):
    """ Find the available viewers for this class. """
    viewers = []
    cls = findClass(className)
    baseClasses = cls.mro()
    for viewer in emViewersDict.values():
        if environment in viewer._environments:
            for t in viewer._targets:
                if t in baseClasses:
                    viewers.append(viewer)
                    break
    return viewers

def findWizards(protocol, environment):
    """ Find availables wizards for this class. 
    Returns:
        a dict with the paramName and wizards for this class."""
    wizards = {}
    baseClasses = protocol.getClass().mro()
    for wiz in emWizardsDict.values():
        if environment in wiz._environments:
            for cls, params in wiz._targets:
                if cls in baseClasses:
                    for p in params:
                        wizards[p] = wiz
    return wizards    

# Update global dictionary with variables found
globals().update(emProtocolsDict)
globals().update(emObjectsDict)
globals().update(emViewersDict)
    