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

from constants import *
from data import *
from data_tiltpairs import *
from protocol import *
from convert import *
from wizards import *
from pyworkflow.utils import importFromPlugin
from pyworkflow.wizard import Wizard
from pyworkflow.viewer import Viewer
from pyworkflow import Config

import pyworkflow.plugin


class Domain(pyworkflow.plugin.Domain):
    _name = __name__
    _objectClass = EMObject
    _protocolClass = Protocol
    _viewerClass = Viewer
    _wizardClass = Wizard
    _baseClasses = globals()


class Plugin(pyworkflow.plugin.Plugin):
    pass


def findClass(className):
    c = Domain.getProtocols().get(
        className,
        Domain.getObjects().get(className, None))

    if c is None:
        raise Exception("findClass: class '%s' not found." % className)

    return c


def findSubClasses(classDict, className):
    """ Find all subclasses of a give className. """
    cls = classDict[className]
    subclasses = {}
    
    for k, v in classDict.iteritems():
        if issubclass(v, cls):
            subclasses[k] = v    
    return subclasses


def getPreferredViewers(className):
    """ Find and import the preferred viewers for this class. """
    preferredViewerNames = Config.VIEWERS.get(className, [])
    if not isinstance(preferredViewerNames, list):
        preferredViewerNames = [preferredViewerNames]
    preferredViewers = []  # we will try to import them and store here
    for prefViewerStr in preferredViewerNames:
        try:
            (prefViewerModule, prefViewerClassName) = prefViewerStr.rsplit('.', 1)
            prefViewer = importFromPlugin(prefViewerModule, prefViewerClassName, doRaise=True)
            preferredViewers.append(prefViewer)
        except Exception as e:
            print("Couldn't load \"%s\" as preferred viewer.\n"
                  "There might be a typo in your VIEWERS "
                  "variable or an error in the viewer's plugin installation" % prefViewerStr)
            print(e)
    return preferredViewers

def findViewers(className, environment):
    """ Find the available viewers for this class. """
    viewers = []
    cls = findClass(className)
    baseClasses = cls.mro()
    preferredViewers = getPreferredViewers(className)
    preferedFlag = 0

    for viewer in Domain.getViewers().values():
        if environment in viewer._environments:
            for t in viewer._targets:
                if t in baseClasses:
                    for prefViewer in preferredViewers:
                        if viewer is prefViewer:
                            viewers.insert(0, viewer)
                            preferedFlag = 1
                            break
                    else:
                        if t == cls:
                            viewers.insert(preferedFlag, viewer)
                        else:
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
    return findWizardsFromDict(protocol, environment, Domain.getWizards())

# Update global dictionary with variables found
#globals().update(emProtocolsDict)
#globals().update(emObjectsDict)
#globals().update(emViewersDict)

def loadSetFromDb(dbName, dbPrefix=''):
    from pyworkflow.mapper.sqlite import SqliteFlatDb
    db = SqliteFlatDb(dbName=dbName, tablePrefix=dbPrefix)
    setClassName = db.getProperty('self') # get the set class name
    setObj = Domain.getObjects()[setClassName](filename=dbName, prefix=dbPrefix)
    return setObj


def runProgram(program, params):
    env = None

    if program.startswith('xmipp'):
        xmipp3 = importFromPlugin('xmipp3', 'Plugin')
        env = xmipp3.getEnviron()
    if program.startswith('relion'):
        relion = importFromPlugin('relion', 'Plugin')
        env = relion.getEnviron()
    elif (program.startswith('e2') or
              program.startswith('sx')):
        eman2 = importFromPlugin('eman2', 'Plugin')
        env = eman2.getEnviron()
    elif program.startswith('b'):
        bsoft = importFromPlugin('bsoft', 'Plugin')
        env = bsoft.getEnviron()

    pwutils.runJob(None, program, params, env=env)
