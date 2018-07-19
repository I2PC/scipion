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
import importlib
import pkgutil
import inspect

import pyworkflow.utils as pwutils
from pyworkflow.utils.reflection import getSubclassesFromModules, getSubclasses
from data import *
from data_tiltpairs import *
from protocol import *
from constants import *
from convert import *
from pyworkflow.wizard import Wizard
from viewer import *
import transformations
import pdb_handler


class PluginMeta(type):
    """ Metaclass that should be used to mark a given Python module
    as a Scipion plugin.
    """
    def __init__(cls, name, bases, dct):
        print("Creating PluginMeta....name: ", name)
        from pyworkflow.utils import prettyDict
        super(PluginMeta, cls).__init__(name, bases, dct)
        # for attrName, attr in cls.__dict__.iteritems():
        #     if '__private' in attrName: # and iscallable(attr):
        #         attr()

        m = Domain.registerPlugin(cls.__module__)
        #
        # if cls.__module__ == 'eman2':
        #     prettyDict(dct)
        #     prettyDict(m.__dict__)


class Domain:
    _plugins = {}
    _protocols = {}
    _objects = {}
    _viewers = {}
    _wizards = {}

    @classmethod
    def __getSubmodule(cls, name, subname):
        try:
            return importlib.import_module('%s.%s' % (name, subname))
        except Exception as e:
            msg = str(e)
            # FIXME: The following is a quick and dirty way to filter
            # when the submodule is not present
            if msg != 'No module named %s' % subname:
                print("  failed to load: %s.%s" % (name, subname))
                print("   error: %s" % e)
            return None

    @classmethod
    def _pluginMethod_hello(cls, plugin):
        print("hello from plugin: ", plugin.__name__)

    @classmethod
    def registerPlugin(cls, name):
        """ Register a new plugin. This function should only be called
        when creating a class with __metaclass__=PluginMeta that will
        trigger this.
        """
        print("Registering plugin: ", name)
        m = importlib.import_module(name)
        cls._plugins[name] = m  # Register the name to as a plugin
        # TODO: Load subclasses (protocols, viewers, wizards)

        # Load bibtex
        m._bibtex = {}
        bib = cls.__getSubmodule(name, 'bibtex')
        if bib is not None:
            try:
                m._bibtex = pwutils.parseBibTex(bib.__doc__)
            except:
                pass

        if not hasattr(m, 'hello'):
            m.hello = lambda: cls._pluginMethod_hello(m)

        return m


    @classmethod
    def __hasPluginMeta(cls, m, className):
        """ Return True if the className refers to a class that
        has PluginMeta as MetaClass (i.e. the module is a plugin)
        """
        a = getattr(m, className)
        return inspect.isclass(a) and type(a) is PluginMeta

    @classmethod
    def __isPlugin(cls, m):
        """ Return True if the module is a Scipion plugin. """
        return m.__name__ in cls._plugins

    @classmethod
    def getPlugins(cls):
        """ Return existing plugins for this Domain. """
        loaded = getattr(cls, '_pluginsLoaded', False)
        if not loaded: # Load plugin only once
            for p, name, isModule in pkgutil.iter_modules():
                if isModule:
                    try:
                        importlib.import_module(name)
                        # NOTE: After importing the modules they will
                        # automatically being register if they have a
                        # class Plugin with metaclass PluginMeta
                    except Exception:
                        pass
            cls._pluginsLoaded = True
        return dict(cls._plugins)

    @classmethod
    def getPlugin(cls, name):
        """ Load a given plugin name. """
        m = importlib.import_module(name)

        if not cls.__isPlugin(m):
            raise Exception("Invalid plugin '%s'. "
                            "Class Plugin with __metaclass__=PluginMeta "
                            "not found" % name)
        return m

    @classmethod
    def __getSubclasses(cls, submoduleName, BaseClass):
        subclasses = getattr(cls, '_%s' % submoduleName)

        if not subclasses:  # Only discover subclasses once
            for pluginName, plugin in cls.getPlugins().iteritems():
                sub = cls.__getSubmodule(pluginName, submoduleName)
                if sub is not None:
                    for name in dir(sub):
                        attr = getattr(sub, name)
                        if inspect.isclass(attr) and issubclass(attr, BaseClass):
                            attr._package = plugin  # Set this special property used by Scipion
                            subclasses[name] = attr
            subclasses.update(getSubclasses(BaseClass, globals()))

        return subclasses

    @classmethod
    def getProtocols(cls):
        """ Return all Protocol subclasses from all plugins for this domain.
        """
        return cls.__getSubclasses('protocols', Protocol)

    @classmethod
    def getObjects(cls):
        """ Return all EMObject subclasses from all plugins for this domain.
        """
        return cls.__getSubclasses('objects', EMObject)

    @classmethod
    def getViewers(cls):
        """ Return all Viewer subclasses from all plugins for this domain.
        """
        return cls.__getSubclasses('viewers', Viewer)

    @classmethod
    def getWizards(cls):
        """ Return all Wizard subclasses from all plugins for this domain.
        """
        return cls.__getSubclasses('wizards', Wizard)

        
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
    env = None

    if program.startswith('xmipp'):
        import xmipp3
        env = xmipp3.getEnviron()
    if program.startswith('relion'):
        import relion
        env = relion.getEnviron()
    elif (program.startswith('e2') or
              program.startswith('sx')):
        import eman2
        env = eman2.getEnviron()
    elif program.startswith('b'):
        import pyworkflow.em.packages.bsoft as bsoft
        env = bsoft.getEnviron()

    pwutils.runJob(None, program, params, env=env)
