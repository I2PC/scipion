# **************************************************************************
# *
# * Authors:     Yaiza Rancel (cyrancel@cnb.csic.es)
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

import sys
import argparse
import os
import re
from future.utils import iteritems

from pyworkflow.plugin import Domain
from plugin_funcs import PluginRepository, PluginInfo
import script

#  ************************************************************************
#  *                                                                      *
#  *                       External (EM) Plugins                          *
#  *                                                                      *
#  ************************************************************************

MODE_INSTALL_PLUGIN = 'installp'
MODE_UNINSTALL_PLUGIN = 'uninstallp'
MODE_LIST_BINS = 'listb'
MODE_INSTALL_BINS = 'installb'
MODE_UNINSTALL_BINS = 'uninstallb'

args = sys.argv[1:]

pluginRepo = PluginRepository()

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
subparsers = parser.add_subparsers(help='mode "installp", "uninstallp" or "listb"',
                                   dest='mode',
                                   title='Mode',
                                   description='available modes are "installp" or "uninstallp"')

############################################################################
#                               Install parser                             #
############################################################################

installParser = subparsers.add_parser("installp", formatter_class=argparse.RawTextHelpFormatter,
                                      usage="%s  [-h] [--noBin] [-p pluginName [pipVersion ...]]" %
                                            (' '.join(args[:2])),
                                      epilog=("Example: %s -p scipion-em-grigoriefflab 1.0.1 "
                                              "-p scipion-em-relion -p scipion-em-eman2 \n\n" %
                                              ' '.join(args[:2])),
                                      add_help=False)
installParser.add_argument('-h', '--help', action='store_true', help='show help')
installParser.add_argument('--noBin', action='store_true',
                            help='Optional flag to install plugins only as a python module,\n'
                                 'without installing the plugin binaries. This will affect\n'
                                 'all plugins specified in the command.')
installParser.add_argument('--checkUpdates', action='store_true',
                           help='Optional flag to check which plugins have new releases.\n')
installParser.add_argument('-p', '--plugin', action='append', nargs='+',
                           metavar=('pluginName', 'pluginVersion'),
                           help='- pluginName:     the name of the plugin to install from the list\n'
                                 '                 of available plugins shown below.\n'
                                 '- pluginVersion: (optional) pip version to install. If not specified,\n'
                                 '                 will install the latest compatible with current Scipion.')

installParser.add_argument('--devel', action='store_true',
                            help='Optional flag to indicate that we will pass install sources instead\n'
                                 'of pip names. Sources might be local paths or git urls. With local\n'
                                 'paths, will do pip install -e (editable mode). It is expected to find\n'
                                 'the plugin name in the basename of the path or in the repo name. \n'
                                 '(i.e. it needs to match the one specified in setup.py). E.g:\n'
                                 'scipion installp -p path/to/pluginName --devel \n'
                                 'scipion installp -p https://github.com/someOrg/pluginName.git --devel')
installParser.add_argument('-j',
                              default='1',
                              metavar='j',
                              help='Number of CPUs to use for compilation \n')

############################################################################
#                             Uninstall parser                             #
############################################################################

uninstallParser = subparsers.add_parser("uninstallp", formatter_class=argparse.RawTextHelpFormatter,
                                        usage="%s  [-h] [-p pluginName [binVersion ...]]" %
                                              (' '.join(args[:2])),
                                        epilog="Example: %s -p scipion-em-grigoriefflab scipion-em-eman2 \n\n" %
                                               ' '.join(args[:2]),
                                        add_help=False)
uninstallParser.add_argument('-h', '--help', action='store_true', help='show help')
uninstallParser.add_argument('--noBin', action='store_true',
                            help='Optional flag to uninstall plugins only as a python module,\n'
                                 'without uninstalling the plugin binaries. This will affect\n'
                                 'all plugins specified in the command.')
uninstallParser.add_argument('-p', '--plugin', action='append',
                             metavar='pluginName',
                             help='The name of the plugin to uninstall from the list\n'
                                  'of available plugins shown below.\n')

############################################################################
#                           Install Bins parser                            #
############################################################################

installBinParser = subparsers.add_parser("installb", formatter_class=argparse.RawTextHelpFormatter,
                                         usage="%s  [-h] binName1 binName2-1.2.3 binName3 ..." %
                                         (' '.join(args[:2])),
                                         epilog="Example: %s ctffind4 unblur-1.0.15\n\n" %
                                         (' '.join(args[:2])),
                                         add_help=False)
#installBinParser.add_argument('pluginName', metavar='pluginName',
#                              help='The name of the plugin whose bins we want to uninstall.\n')
installBinParser.add_argument('-h', '--help', action='store_true', help='show help')
installBinParser.add_argument('binName', nargs='*',
                              metavar='binName(s)',
                              help='The name(s) of the bins we want install, optionally with \n'
                                   'version in the form name-version. If no version is specified,\n'
                                   'will install the last one.')
installBinParser.add_argument('-j',
                              default='1',
                              metavar='j',
                              help='Number of CPUs to use for compilation \n')

############################################################################
#                          Uninstall Bins parser                           #
############################################################################

uninstallBinParser = subparsers.add_parser("uninstallb", formatter_class=argparse.RawTextHelpFormatter,
                                           usage="%s [-h] binName1 binName2-1.2.3 binName3 ..." %
                                           (' '.join(args[:2])),
                                           epilog="Example: %s ctffind4 unblur-1.0.15\n\n " %
                                           (' '.join(args[:2])),
                                           add_help=False)
#uninstallBinParser.add_argument('pluginName', metavar='pluginName',
#                                help='The name of the plugin whose bins we want to uninstall.\n')
uninstallBinParser.add_argument('-h', '--help', action='store_true', help='show help')
uninstallBinParser.add_argument('binName', nargs='+',
                                metavar='binName(s)',
                                help='The name(s) of the bins we want to uninstall\n'
                                     '(optionally with version in the form name-version). \n'
                                     'If no version is specified, will uninstall the last one.\n')

modeToParser = {MODE_INSTALL_BINS: installBinParser,
                MODE_UNINSTALL_BINS: uninstallBinParser,
                MODE_INSTALL_PLUGIN: installParser,
                MODE_UNINSTALL_PLUGIN: uninstallParser}

parsedArgs = parser.parse_args(args[1:])
mode = parsedArgs.mode
parserUsed = modeToParser[mode]
exitWithErrors = False

if parsedArgs.help or (mode in [MODE_INSTALL_BINS, MODE_UNINSTALL_BINS]
                       and len(parsedArgs.binName) == 0):

    if mode not in [MODE_INSTALL_BINS, MODE_UNINSTALL_BINS]:
        parserUsed.epilog += pluginRepo.printPluginInfoStr()
    else:
        env = script.defineBinaries([])
        env.setDefault(False)
        installedPlugins = Domain.getPlugins()
        for p, pobj in iteritems(installedPlugins):
            pobj.Plugin.defineBinaries(env)
        parserUsed.epilog += env.printHelp()
    parserUsed.print_help()
    parserUsed.exit(0)

elif mode == MODE_INSTALL_PLUGIN:
    if parsedArgs.checkUpdates:
        print(pluginRepo.printPluginInfoStr(withUpdates=True))
        installParser.exit(0)

    if parsedArgs.devel:
        for p in parsedArgs.plugin:
            pluginSrc = p[0]
            pluginName = ""
            if os.path.exists(pluginSrc):
                pluginName = os.path.basename(pluginSrc.rstrip('/'))
                numberProcessor = parsedArgs.j
            else:  # we assume it is a git url
                m = re.match('https://github.com/(.*)/(.*).git', pluginSrc)
                if m:
                    pluginName = m.group(2)
            if not pluginName:
                print("ERROR: Couldn't find pluginName for source %s" % pluginSrc)
                exitWithErrors = True
            else:
                plugin = PluginInfo(pipName=pluginName, pluginSourceUrl=pluginSrc, remote=False)
                processors = parsedArgs.j
                installed = plugin.installPipModule()
                if installed and not parsedArgs.noBin:
                    plugin.installBin(args=['-j', numberProcessor])
    else:
        pluginsToInstall = list(zip(*parsedArgs.plugin))[0]
        pluginDict = pluginRepo.getPlugins(pluginList=pluginsToInstall,
                                           getPipData=True)
        if not pluginDict:
            exitWithErrors = True
        else:
            for cmdTarget in parsedArgs.plugin:
                pluginName = cmdTarget[0]
                pluginVersion = "" if len(cmdTarget) == 1 else cmdTarget[1]
                numberProcessor = parsedArgs.j
                plugin = pluginDict.get(pluginName, None)
                if plugin:
                    installed = plugin.installPipModule(version=pluginVersion)
                    if installed and not parsedArgs.noBin:
                        plugin.installBin(args=['-j', numberProcessor])
                else:
                    print("WARNING: Plugin %s does not exist." % pluginName)
                    exitWithErrors = True

elif parsedArgs.mode == MODE_UNINSTALL_PLUGIN:

    for pluginName in parsedArgs.plugin:
        plugin = PluginInfo(pluginName, pluginName, remote=False)
        if plugin.isInstalled():
            if not parsedArgs.noBin:
                plugin.uninstallBins()
            plugin.uninstallPip()
        else:
            print("WARNING: Plugin %s is not installed." % pluginName)

elif parsedArgs.mode == MODE_INSTALL_BINS:

    binToInstallList = parsedArgs.binName
    binToPlugin = pluginRepo.getBinToPluginDict()
    for binTarget in binToInstallList:
        pluginTargetName = binToPlugin.get(binTarget, None)
        if pluginTargetName is None:
            print('ERROR: Could not find target %s' % binTarget)
            continue
        pmodule = Domain.getPlugin(pluginTargetName)
        numberProcessor = parsedArgs.j
        pinfo = PluginInfo(name=pluginTargetName, plugin=pmodule, remote=False)
        pinfo.installBin([binTarget, '-j', numberProcessor])


elif parsedArgs.mode == MODE_UNINSTALL_BINS:

    binToInstallList = parsedArgs.binName
    binToPlugin = pluginRepo.getBinToPluginDict()
    for binTarget in binToInstallList:
        pluginTargetName = binToPlugin.get(binTarget, None)
        if pluginTargetName is None:
            print('ERROR: Could not find target %s' % binTarget)
            continue
        pmodule = Domain.getPlugin(pluginTargetName)
        pinfo = PluginInfo(name=pluginTargetName, plugin=pmodule, remote=False)
        pinfo.uninstallBins([binTarget])

if exitWithErrors:
    parserUsed.exit(1)
else:
    parserUsed.exit(0)
