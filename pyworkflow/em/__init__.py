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
This modules contains classes related with EM
"""
#import os
#import sys

#import pyworkflow as pw
#import pyworkflow.utils as pwutils
from pyworkflow.utils.reflection import getSubclassesFromModules, getSubclasses, getModules
from data import *
from data_tiltpairs import *
from protocol import *
from constants import *
from convert import *
from pyworkflow.wizard import Wizard
from viewer import *
import transformations
import pdb_handler

PACKAGES_PATH = os.path.join(pw.HOME, 'em', 'packages')
_emPackagesDict = None

def getPackages():
    global _emPackagesDict
    if _emPackagesDict is None:
        sys.path.insert(0, PACKAGES_PATH)
        _emPackagesDict = getModules(PACKAGES_PATH)
        sys.path.pop(0)
    return _emPackagesDict

# Load all Protocol subclasses found in EM-packages
_emProtocolsDict = None

def getProtocols():
    """ Load all protocols subclasses defined in all em-packages. """
    global _emProtocolsDict
    if _emProtocolsDict is None:
        _emProtocolsDict = getSubclassesFromModules(Protocol, getPackages())
        _emProtocolsDict.update(getSubclasses(Protocol, globals()))
    return _emProtocolsDict

_emObjectsDict = None 

def getObjects():
    """ Load all EMObject subclasses found in EM-packages. """
    global _emObjectsDict
    if _emObjectsDict is None:        
        _emObjectsDict = getSubclassesFromModules(EMObject, getPackages())
        _emObjectsDict.update(getSubclasses(EMObject, globals()))
    return _emObjectsDict

_emViewersDict = None

def getViewers():
    """ Load all subclasses of Viewer of different packages. """
    global _emViewersDict
    if _emViewersDict is None:
        _emViewersDict = getSubclassesFromModules(Viewer, getPackages())
    return _emViewersDict

_emWizardsDict = None

def getWizards():
    """ Load all subclasses of Wizards. """
    global _emWizardsDict
    if _emWizardsDict is None:
        _emWizardsDict = getSubclassesFromModules(Wizard, getPackages())
    return _emWizardsDict
        
        
def findClass(className):
    
    if className in getProtocols():
        return getProtocols()[className]
    
    if className in getObjects():
        return getObjects()[className]
    
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
    for viewer in getViewers().values():
        if environment in viewer._environments:
            for t in viewer._targets:
                if t in baseClasses:
                    viewers.append(viewer)
                    break
    return viewers


#TODO: If we divide the way to find wizards for web
# and desktop, the environment is no longer needed
def findWizardsFromDict(protocol, environment, wizDict):
    wizards = {}
    baseClasses = [cls.__name__ for cls in protocol.getClass().mro()]
    
    for wiz in wizDict.values():
        if environment in wiz._environments:
            for cls, params in wiz._targets:
                if cls.__name__ in baseClasses:
                    for p in params:
                        wizards[p] = wiz
    return wizards

    
def findWizards(protocol, environment):
    """ Find availables wizards for this class. 
    Returns:
        a dict with the paramName and wizards for this class."""
    return findWizardsFromDict(protocol, environment, getWizards())

# Update global dictionary with variables found
#globals().update(emProtocolsDict)
#globals().update(emObjectsDict)
#globals().update(emViewersDict)

def loadSetFromDb(dbName, dbPrefix=''):
    from pyworkflow.mapper.sqlite import SqliteFlatDb
    db = SqliteFlatDb(dbName=dbName, tablePrefix=dbPrefix)
    setClassName = db.getProperty('self') # get the set class name
    setObj = getObjects()[setClassName](filename=dbName, prefix=dbPrefix)
    return setObj


def runProgram(program, params):
    env = None

    if program.startswith('xmipp'):
        import pyworkflow.em.packages.xmipp3 as xmipp3
        env = xmipp3.getEnviron()
    if program.startswith('relion'):
        import pyworkflow.em.packages.relion as relion
        env = relion.getEnviron()
    elif (program.startswith('e2') or
              program.startswith('sx')):
        import pyworkflow.em.packages.eman2 as eman2
        env = eman2.getEnviron()
    elif program.startswith('b'):
        import pyworkflow.em.packages.bsoft as bsoft
        env = bsoft.getEnviron()

    pwutils.runJob(None, program, params, env=env)