#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] MRC Laboratory of Molecular Biology (MRC-LMB)
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

from __future__ import print_function
import sys
import importlib
import inspect
import traceback
from collections import OrderedDict
import argparse

from pyworkflow.em import Domain
from pyworkflow.protocol import Protocol
from install.plugin_funcs import PluginInfo
import pyworkflow.em as em
import pyworkflow.utils as pwutils
from install.funcs import Environment


def getSubmodule(name, subname):
    """ Return a tuple: (module, error)
    If module is None:
        1) if error is None is that the submodule does not exist
        2) if error is not None, it is an Exception raised when
        importing the submodule
    """
    try:
        m = importlib.import_module('%s.%s' % (name, subname))
        r = (m, None)
    except Exception as e:
        noModuleMsg = 'No module named %s' % subname
        msg = str(e)
        r = (None, None if msg == noModuleMsg else traceback.format_exc())
    return r


def getFirstLine(doc):
    """ Get the first non empty line from doc. """
    if doc:
        for lines in doc.split('\n'):
            l = lines.strip()
            if l:
                return l
    return ''


def checkPlugins():
    """ Discover all plugins and print some info.
    """
    plugins = Domain.getPlugins()
    env = Environment()
    # pip and python are needed by Xmipp plugin
    env.addTarget('pip')
    env.addTarget('python')

    print("Plugins:")
    for k, plugin in plugins.iteritems():
        print("-", plugin.__name__)
        plugin.Plugin.defineBinaries(env)

    print("\nDefined binaries: ")
    print(env.printHelp())

    if args.verbosity > 0:
        print("Objects")
        pwutils.prettyDict(Domain.getObjects())

        print("Protocols")
        pwutils.prettyDict(Domain.getProtocols())

        print("Viewers")
        pwutils.prettyDict(Domain.getViewers())


def checkPlugin(pluginName):
    """ Check submodules of a given plugin. """
    for subName in ['constants', 'convert', 'protocols',
                    'wizards', 'viewers', 'tests']:
        sub, error = getSubmodule(pluginName, subName)

        if sub is None:
            if error is None:
                msg = " missing"
            else:
                msg = " error -> %s" % error

        else:
            msg = " loaded"

        print("   >>> %s: %s" % (subName, msg))


def printPluginInfo(pluginName, verbosity):
    """ Print info about a given plugin """
    showBase = verbosity > 0
    subclasses = {}
    emCategories = [('Imports', em.ProtImport),
                    ('Micrographs', em.ProtMicrographs),
                    ('Particles', em.ProtParticles),
                    ('2D', em.Prot2D),
                    ('3D', em.Prot3D)]

    plugin = Domain.getPlugin(pluginName)
    version = PluginInfo('scipion-em-%s' % pluginName).pipVersion
    print("Plugin name: %s, version: %s\n" % (pluginName, version))

    print("Plugin binaries: ")
    env = Environment()
    plugin.Plugin.defineBinaries(env)
    print(env.printHelp())

    # print bibtex
    bib, error2 = getSubmodule(pluginName, 'bibtex')
    if bib is None:
        if error2 is None:
            msg = " missing bibtex"
        else:
            msg = " error -> %s" % error2
    else:
        print("Plugin references:")
        bibtex = pwutils.parseBibTex(bib.__doc__)

        for citeStr in bibtex:
            text = Protocol()._getCiteText(bibtex[citeStr])
            print(text)

    # print protocols
    sub, error = getSubmodule(pluginName, 'protocols')
    if sub is None:
        if error is None:
            msg = " missing protocols"
        else:
            msg = " error -> %s" % error

    else:
        for name in dir(sub):
            attr = getattr(sub, name)
            if inspect.isclass(attr) and issubclass(attr, Protocol):
                # Set this special property used by Scipion
                attr._package = plugin
                subclasses[name] = attr

    print("Plugin protocols:\n")
    print("%-35s %-35s %-10s %-s" % (
        'NAME', 'LABEL', 'CATEGORY', 'DESCRIPTION'))

    prots = OrderedDict(sorted(subclasses.items()))
    for prot in prots:
        label = prots[prot].getClassLabel()
        desc = getFirstLine(prots[prot].__doc__)
        cat = 'None'

        for c in emCategories:
            if issubclass(prots[prot], c[1]):
                cat = c[0]
            if prots[prot].isBase():
                cat = 'Base prot'

        # skip Base protocols if not requested
        if prots[prot].isBase() and not showBase:
            continue
        else:
            print("%-35s %-35s %-10s %-s" % (prot, label, cat, desc))


def installBinary(plugin, target):
    env = Environment()
    plugin.Plugin.defineBinaries(env)

    if env.hasTarget(target):
        env.executeTargets(target)
    else:
        print("ERROR: Target %s not found, please provide a valid"
              "version from the ones listed below" % target)
        print(env.printHelp())


parser = argparse.ArgumentParser(description="Script to provide information about installed plugins.")
add = parser.add_argument  # shortcut
add("plugin", help="Plugin name", nargs='?', default="")
add("--check", action="store_true",
    help="Check plugin code and submodules. "
         "If verbosity > 0, all discovered objects, protocols and viewers will be shown. ")
add("--info", action="store_true",
    help="Show plugin summary (protocols, citations, etc). "
         "If verbosity > 0, base protocol classes will be shown. ")
add("--install", nargs='+',
    help="Specify a binary version to install. "
         "If you pass a second argument, it should a path to the installation"
         "and a link will be created. ")
add("-v", "--verbosity", action="count", default=0)
args = parser.parse_args()

pluginName = args.plugin

if not pluginName:  # List all plugins
    checkPlugins()
else:
    plugin = Domain.getPlugin(pluginName)
    print("Plugin: %s" % pluginName)

    if args.install:
        if len(args.install) > 2:
            print("ERROR: --install option should have one or two arguments. ")
        installBinary(plugin, target=args.install[0])
        if len(args.install) > 1:
            print("TODO: Installation via a symbolic link is not yet implemented. ")
    else:
        if args.check:
            checkPlugin(pluginName)

        if args.info:
            printPluginInfo(pluginName, args.verbosity)
