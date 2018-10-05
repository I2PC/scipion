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
from install.plugin_funcs import PluginRepository, PluginInfo
import argparse

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
                                      epilog="Example: %s -p grigoriefflab 1.0.1 -p relion -p eman \n\n" %
                                      ' '.join(args[:2]),
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

############################################################################
#                             Uninstall parser                             #
############################################################################

uninstallParser = subparsers.add_parser("uninstallp", formatter_class=argparse.RawTextHelpFormatter,
                                        usage="%s  [-h] [-p pluginName [binVersion ...]]" %
                                              (' '.join(args[:2])),
                                        epilog="Example: %s grigoriefflab eman \n\n" %
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
                                         epilog="Example: %s ctffind4 unblur-1.0.15\n\n %s" %
                                             (' '.join(args[:2]), pluginRepo.printPluginInfoStr(withBins=True)))
#installBinParser.add_argument('pluginName', metavar='pluginName',
#                              help='The name of the plugin whose bins we want to uninstall.\n')
installBinParser.add_argument('binName', nargs='+',
                              metavar='binName(s)',
                              help='The name(s) of the bins we want install, optionally with \n'
                                   'version in the form name-version. If no version is specified,\n'
                                   'will install the last one.')
installBinParser.add_argument('-j',
                              default=1,
                              metavar='j',
                              help='Number of CPUs to use for compilation \n')

############################################################################
#                          Uninstall Bins parser                           #
############################################################################

uninstallBinParser = subparsers.add_parser("uninstallb", formatter_class=argparse.RawTextHelpFormatter,
                                           usage="%s [-h] binName1 binName2-1.2.3 binName3 ..." %
                                           (' '.join(args[:2])),
                                           epilog="Example: %s ctffind4 unblur-1.0.15\n\n %s" %
                                           (' '.join(args[:2]), pluginRepo.printPluginInfoStr(withBins=True)))
#uninstallBinParser.add_argument('pluginName', metavar='pluginName',
#                                help='The name of the plugin whose bins we want to uninstall.\n')
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

if mode not in [MODE_INSTALL_BINS, MODE_UNINSTALL_BINS] and parsedArgs.help:
    parserUsed = modeToParser[mode]
    parserUsed.epilog += pluginRepo.printPluginInfoStr()
    parserUsed.print_help()
    parserUsed.exit(0)

elif mode == MODE_INSTALL_PLUGIN:
    if parsedArgs.checkUpdates:
        print(pluginRepo.printPluginInfoStr(withUpdates=True))
        installParser.exit(0)

    pluginDict = pluginRepo.getPlugins(pluginList=list(zip(*parsedArgs.plugin))[0],
                                       getPipData=True)
    if not pluginDict:
        print('\n' + installParser.epilog)
    else:
        for cmdTarget in parsedArgs.plugin:
            pluginName = cmdTarget[0]
            pluginVersion = "" if len(cmdTarget) == 1 else cmdTarget[1]
            plugin = pluginDict.get(pluginName, None)
            if plugin:
                installed = plugin.installPipModule(version=pluginVersion)
                if installed and not parsedArgs.noBin:
                    plugin.installBin()
            else:
                print("WARNING: Plugin %s does not exist." % pluginName)

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

    # Get all plugins installed
    pluginDict = pluginRepo.getPlugins(
        pluginList=list(zip(*parsedArgs.binName))[0],
        getPipData=True)
    binToInstallList = parsedArgs.binName

    if not pluginDict:
        print('\n' + installParser.epilog)
    else:
        # Move through the plugin list
        for pluginName in pluginDict:
            if binToInstallList:
                plugin = pluginDict.get(pluginName, None)
                if plugin.isInstalled():
                    # Get plugin binaries
                    binariesPluginList = plugin.getBinVersions()
                    # For each binary to install
                    for binName in binToInstallList:
                        binVersions = [binary for binary in binariesPluginList
                                       if binary.split('-')[0] == binName or
                                       binary == binName]
                        if binVersions:
                            if binName not in binVersions:
                                binToInstallName = binVersions[len(binVersions)-1]
                            else:
                                binToInstallName = binName
                            try:
                                plugin.installBin(args=[binToInstallName])
                                binToInstallList.remove(binName)
                            except AssertionError as err:  # TODO The correct exception must be captured
                                print("WARNING: Binaries of %s has not been "
                                      "installed." % binName)
                                print("WARNING: Binaries of %s does not exist."
                                      % binName)
                                binToInstallList.remove(binName)
            else:
                break
        if binToInstallList:
            print("\n----------------------------------- \n")
            for binName in binToInstallList:
                print("WARNING: Binaries of %s has not been installed."
                      % binName)
                print("WARNING: Binaries of %s does not exist. \n"
                      % binName)
            print("----------------------------------- \n")

elif parsedArgs.mode == MODE_UNINSTALL_BINS:

    # Get all plugins installed
    pluginDict = pluginRepo.getPlugins(
        pluginList=list(zip(*parsedArgs.binName))[0],
        getPipData=True)
    binToUninstallList = parsedArgs.binName

    if not pluginDict:
        print('\n' + installParser.epilog)
    else:
        # Move through the plugin list
        for pluginName in pluginDict:
            if binToUninstallList:
                plugin = pluginDict.get(pluginName, None)
                if plugin.isInstalled():
                    # Get plugin binaries
                    binariesPluginList = plugin.getBinVersions()
                    # For each binary to install
                    for binName in binToUninstallList:
                        binVersions = [binary for binary in binariesPluginList
                                       if binary.split('-')[0] == binName or
                                       binary == binName]
                        if binVersions:
                            if binName not in binVersions:
                                binToInstallName = [binVersions[len(binVersions)-1]]
                            else:
                                binToInstallName = [binName]
                            try:
                                plugin.uninstallBins(binToInstallName)
                                print("Binaries of %s have been uninstalled "
                                      "successfully." % binName)
                                binToUninstallList.remove(binName)
                                
                            except AssertionError as err:  # TODO The correct exception must be captured
                                print("WARNING: Binaries of %s have not been "
                                      "uninstalled." % binName)
                                print("WARNING: Binaries of %s don't exist."
                                      % binName)
                                binToUninstallList.remove(binName)
            else:
                break
        if binToUninstallList:
            print("\n----------------------------------- \n")
            for binName in binToUninstallList:
                print("WARNING: Binaries of %s have not been uninstalled."% binName)
                print("WARNING: Binaries of %s don't exist. \n" % binName)
            print("----------------------------------- \n")
