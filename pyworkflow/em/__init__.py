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
from pyworkflow.utils.reflection import getSubClassesFromPath, getSubclasses
from data import *
from protocol import *
from viewer import Viewer
#from packages import *

PACKAGES_PATH = os.path.join(pw.HOME, 'em', 'packages')

# Load all Protocol subclasses found in EM-packages
emProtocolsDict = getSubClassesFromPath(Protocol, PACKAGES_PATH)
emProtocolsDict.update(getSubclasses(Protocol, globals()))

# Load all EMObject subclasses found in EM-packages
emObjectsDict = getSubClassesFromPath(EMObject, PACKAGES_PATH)
emObjectsDict.update(getSubclasses(EMObject, globals()))

# Load all subclasses of Viewer of different packages
emViewersDict = getSubClassesFromPath(Viewer, PACKAGES_PATH)
## Get for which objects are viewers registered
#emViewerTargets = {}
#
#for viewer in emViewersDict.values():
#    print "viewer: ", viewer
#    for t in viewer._targets:
#        print "   target: ", t
#        if not t in emViewerTargets:
#            emViewerTargets[t] = []
#        emViewerTargets[t].append(viewer)
        
        
def findClass(className):
    if className in emProtocolsDict:
        return emProtocolsDict[className]
    if className in emObjectsDict:
        return emObjectsDict[className]
    raise Exception("findClass: class '%s' not found." % className)
    
def findViewers(className, environment):
    """ Find the available viewers for this class. """
    viewers = []
    cls = findClass(className)
    baseClasses = cls.mro()
    for viewer in emViewersDict.values():
        if viewer._environment == environment:
            for t in viewer._targets:
                if t in baseClasses:
                    viewers.append(viewer)
                    break
    return viewers

# Update global dictionary with variables found
globals().update(emProtocolsDict)
globals().update(emObjectsDict)
globals().update(emViewersDict)
    