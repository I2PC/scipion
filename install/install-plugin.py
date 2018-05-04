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

MODE_INSTALL_PLUGIN = 'install_plugin'
MODE_UNINSTALL_PLUGIN = 'uninstall_plugin'

args = sys.argv[1:]

pluginRepo = PluginRepository()

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
subparsers = parser.add_subparsers(help='mode "install_plugin" or "uninstall_plugin"',
                                   dest='mode',
                                   title='Mode',
                                   description='available modes are "install_plugin" or "uninstall_plugin"')

installParser = subparsers.add_parser("install_plugin", formatter_class=argparse.RawTextHelpFormatter,
                                      usage="%s  [-h] [--noBin] [-p pluginName [binVersion ...]]" %
                                            (' '.join(args[:2])),
                                      epilog="Example: %s -p ctffind 2.0.4 -p relion -p eman \n\n" %
                                      ' '.join(args[:2]),
                                      add_help=False)
installParser.add_argument('-h', '--help', action='store_true', help='show help')


installParser.add_argument('--noBin', action='store_true',
                            help='Optional flag to install plugins only as a python module,\n'
                                 'without installing the plugin binaries. This will affect\n'
                                 'all plugins specified in the command.')

installParser.add_argument('-p', '--plugin', action='append', nargs='+',
                           metavar=('pluginName', 'pluginVersion'),
                           help='- pluginName: the name of the plugin to install from the list\n'
                                 '             of available plugins shown below.\n'
                                 '- pluginVersion: (optional) pip version to install')

uninstallParser = subparsers.add_parser("uninstall_plugin", formatter_class=argparse.RawTextHelpFormatter,
                                        usage="%s  [-h] [-p pluginName [binVersion ...]]" %
                                              (' '.join(args[:2])),
                                        epilog="Example: %s ctffind eman \n\n" %
                                               ' '.join(args[:2]),
                                        add_help=False)
uninstallParser.add_argument('-h', '--help', action='store_true', help='show help')
uninstallParser.add_argument('--noBin', action='store_true',
                            help='Optional flag to uninstall plugins only as a python module,\n'
                                 'without uninstalling the plugin binaries. This will affect\n'
                                 'all plugins specified in the command.')
uninstallParser.add_argument('-p', '--plugin', action='append',
                             metavar='pluginName',
                             help='- pluginName: the name of the plugin to uninstall from the list\n'
                                  '             of available plugins shown below.\n')


parsedArgs = parser.parse_args(args[1:])
mode = parsedArgs.mode
if parsedArgs.help:
    parserUsed = installParser if mode == MODE_INSTALL_PLUGIN else uninstallParser
    parserUsed.epilog += pluginRepo.printPluginInfo()
    parserUsed.print_help()
    parserUsed.exit(0)
elif mode == MODE_INSTALL_PLUGIN:
    pluginDict = pluginRepo.getPlugins(pluginList=list(zip(*parsedArgs.plugin))[0], getPipData=True)
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

elif parsedArgs.mode == MODE_UNINSTALL_PLUGIN:
    if parsedArgs.help:
        uninstallParser.epilog += pluginRepo.printPluginInfo()
        uninstallParser.print_help()
        uninstallParser.exit(0)
    for pluginName in parsedArgs.plugin:
        plugin = PluginInfo(pluginName, remote=False)
        if plugin.isInstalled():
            if not parsedArgs.noBin:
                plugin.uninstallBins()
            plugin.uninstallPip()
        else:
            print("Plugin %s is not installed" % pluginName)
