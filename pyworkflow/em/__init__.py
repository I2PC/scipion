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
from pyworkflow.utils import importFromPlugin
from pyworkflow.wizard import Wizard
from pyworkflow.viewer import Viewer

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


def findViewers(className, environment):
    """ Find the available viewers for this class. """
    viewers = []
    cls = findClass(className)
    baseClasses = cls.mro()
    for viewer in Domain.getViewers().values():
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
    plugin = None

    if program.startswith('xmipp'):
        plugin = importFromPlugin('xmipp3', 'Plugin')
    elif program.startswith('relion'):
        plugin = importFromPlugin('relion', 'Plugin')
    elif program.startswith('e2') or program.startswith('sx'):
        plugin = importFromPlugin('eman2', 'Plugin')
    elif program.startswith('b'):
        plugin = importFromPlugin('bsoft', 'Plugin')

    env = plugin.getEnviron() if plugin is not None else None

    pwutils.runJob(None, program, params, env=env)
