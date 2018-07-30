# **************************************************************************
# *
# * Authors:     Yaiza Rancel (cyrancel@cnb.csic.es) [1]
# *              Pablo Conesa (pconse@cnb.csic.es) [1]
# *              J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] SciLifeLab, Stockholm University
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

import os
import importlib
import pkgutil
import inspect
import pkg_resources
from email import message_from_string
from collections import OrderedDict
from abc import ABCMeta, abstractmethod

import pyworkflow.utils as pwutils


class Domain:
    """
    Class to represent the application domain.
    It will allow to specify new objects, protocols, viewers and wizards
    through the registration of new plugins.
    """

    # The following classes should be defined in subclasses of Domain
    _protocolClass = None
    _objectClass = None
    _viewerClass = None
    _wizardClass = None
    _baseClasses = {}  # Update this with the base classes of the Domain

    # Dictionaries to store different type of objects
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
            except Exception:
                pass

        return m

    @classmethod
    def __isPlugin(cls, m):
        """ Return True if the module is a Scipion plugin. """
        return m.__name__ in cls._plugins

    @classmethod
    def getPlugins(cls):
        """ Return existing plugins for this Domain. """
        loaded = getattr(cls, '_pluginsLoaded', False)
        if not loaded:  # Load plugin only once
            for p, name, isModule in pkgutil.iter_modules():
                if isModule:
                    try:
                        importlib.import_module(name)
                        # NOTE: After importing the modules they will
                        # automatically being register if they has
                        # called Domain.registerPlugin(__name__)
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
                            # Set this special property used by Scipion
                            attr._package = plugin
                            subclasses[name] = attr
            subclasses.update(
                pwutils.getSubclasses(BaseClass, cls._baseClasses))

        return subclasses

    @classmethod
    def getProtocols(cls):
        """ Return all Protocol subclasses from all plugins for this domain.
        """
        return cls.__getSubclasses('protocols', cls._protocolClass)

    @classmethod
    def getObjects(cls):
        """ Return all EMObject subclasses from all plugins for this domain.
        """
        return cls.__getSubclasses('objects', cls._objectClass)

    @classmethod
    def getViewers(cls):
        """ Return all Viewer subclasses from all plugins for this domain.
        """
        return cls.__getSubclasses('viewers', cls._viewerClass)

    @classmethod
    def getWizards(cls):
        """ Return all Wizard subclasses from all plugins for this domain.
        """
        return cls.__getSubclasses('wizards', cls._wizardClass)


class Plugin:
    __metaclass__ = ABCMeta

    _vars = {}
    _homeVar = ''  # Change in subclasses to define the "home" variable
    _pathVars = []
    _supportedVersions = []

    @classmethod
    def _defineVar(cls, varName, defaultValue):
        """ Internal method to define variables from the environment. """
        cls._vars[varName] = os.environ.get(varName, defaultValue)

    @classmethod
    def _defineEmVar(cls, varName, defaultValue):
        """ Shortcut method to define variables by prepending EM_ROOT
        to the default value.
        """
        cls._defineVar(varName,
                       os.path.join(os.environ['EM_ROOT'], defaultValue))

    @classmethod
    @abstractmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch programs. """
        pass

    @classmethod
    @abstractmethod
    def _defineVariables(cls):
        """ Method to define variables and their default values.
        It will use the method _defineVar that will take a variable value
        from the environment or from an optional default value.

        This method is not supposed to be called from client code,
        except from the Domain class when registering a Plugin.
        """
        pass

    @classmethod
    @abstractmethod
    def defineBinaries(cls, env):
        """ Define required binaries in the given Environment. """
        pass

    @classmethod
    def getVar(cls, varName, defaultValue=None):
        """ Return the value of a given variable. """
        return cls._vars.get(varName, defaultValue)

    @classmethod
    def getHome(cls, *paths):
        """ Return a path from the "home" of the package
         if the _homeVar is defined in the plugin. """
        home = cls.getVar(cls._homeVar)
        return os.path.join(home, *paths) if home else ''

    @classmethod
    def getSupportedVersions(cls):
        """ Return the list of supported binary versions. """
        return cls._supportedVersions

    @classmethod
    def getActiveVersion(cls):
        """ Return the version of the Relion binaries that is currently active.
        In the current implementation it will be inferred from the RELION_HOME
        variable, so it should contain the version number in it. """
        home = cls.getHome()
        for v in cls.getSupportedVersions():
            if v in home:
                return v
        return ''

    @classmethod
    def validateInstallation(cls):
        """
        Check if the binaries are properly installed and if not, return
        a list with the error messages.

        The default implementation will check if the _pathVars exists.
        """
        environ = cls.getEnviron()
        missing = ["%s: %s" % (var, environ[var])
                   for var in cls._pathVars if not os.path.exists(environ[var])]

        return (["Missing variables:"] + missing) if missing else []


class PluginInfo:
    """
    Information related to a given plugin when it is installed via PIP
    """
    def __init__(self, name):
        try:
            dist = pkg_resources.get_distribution(name)
            lines = [l for l in dist._get_metadata(dist.PKG_INFO)]
            tuples = message_from_string('\n'.join(lines))

        except Exception:
            print("Plugin %s seems is not a pip module yet. "
                  "No metadata found" % name)
            tuples = message_from_string('Author: plugin in development mode?')

        self._name = name
        self._metadata = OrderedDict()

        for v in tuples.items():
            if v[0] == 'Keywords':
                break
            self._metadata[v[0]] = v[1]

    def getAuthor(self):
        return self._metadata.get('Author', "")

    def getAuthorEmail(self):
        return self._metadata.get('Author-email', '')

    def getHomePage(self):
        return self._metadata.get('Home-page', '')

    def getKeywords(self):
        return self._metadata.get('Keywords', '')


# FIXME: Remove the following OLD plugin class when not longer needed
# class Plugin(object):
#
#     def __init__(self, name, version=None, configVars=None,
#                  logo=None, bibtex=None):
#         # Plugin name - only mandatory attribute
#         self.name = name
#         # Plugin version
#         self.version = version
#         # dict with default values for env vars needed
#         self.configVars = configVars
#         # path to the logo, relative to each plugin's plugin.py file
#         self.logo = logo
#         # List with the default plugin references e.g. []
#         self.bibtex = bibtex
#         # Set default env vars
#         self.setDefaultEnviron()
#         # plugin pip metadata
#         self.metadata = self.getMetadata()
#
#     def setDefaultEnviron(self):
#         for k in self.configVars:
#             os.environ.setdefault(k, self.configVars[k])
#         environ = pwutils.Environ(os.environ)
#         return environ
#
#     def registerPluginBinaries(self, env):
#         """Overwrite in subclass"""
#         pass
#
#     def getMetadata(self):
#         try:
#             pipPackage = pkg_resources.get_distribution(self.name)
#             metadataLines = [l for l in pipPackage._get_metadata(pipPackage.PKG_INFO)]
#             metadataTuples = message_from_string('\n'.join(metadataLines))
#
#         except Exception as e:
#             print("Plugin %s seems is not a pip module yet. No metadata found" % self.name)
#             metadataTuples = message_from_string('Author: plugin in development mode?')
#
#         metadata = OrderedDict()
#         for v in metadataTuples.items():
#             if v[0] == 'Keywords':
#                 break
#             metadata[v[0]] = v[1]
#         return metadata
#
#     def getAuthor(self):
#         return self.metadata.get('Author', "")
#
#     def getAuthorEmail(self):
#         return self.metadata.get('Author-email', '')
#
#     def getHomePage(self):
#         return self.metadata.get('Home-page', '')
#
#     def getKeywords(self):
#         return self.metadata.get('Keywords', '')
