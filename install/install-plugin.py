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
from install.funcs import Environment
from install.plugin_funcs import PluginRepository
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
                                      epilog="Example: %s -p ctffind -p relion 0 -p eman 2.11  \n\n%s" %
                                      (' '.join(args[:2]), pluginRepo.printPluginInfo()))
installParser.add_argument('--noBin', action='store_true',
                            help='Optional flag to install plugins only as a python module,\n'
                                 'without installing the plugin binaries. This will affect\n'
                                 'all plugins specified in the command.')
installParser.add_argument('-p', '--plugin', action='append', nargs='+',
                           metavar=('pluginName', 'binVersion'),
                           help='- pluginName: the name of the plugin to install from the list\n'
                                 '              of available plugins shown below.\n'
                                 '- binVersion: is an optional param to control the binary version\n'
                                 '              installed. Its possible values are:\n'
                                 "              - None: latest version will be automatically installed\n"
                                 "                      (See example for ctffind below).\n"
                                 "              - 0:    prevents the installation of the plugin binaries\n"
                                 "                      (See example for relion below)\n"
                                 "              - N:    will install the specified version. Available\n"
                                 "                      (See example for eman below)\n"
                                 "Several plugins may be specified in a single command, \n"
                                 'each following the "-p pluginName [binVersion]" pattern. \n')

uninstallParser = subparsers.add_parser("uninstall_plugin", formatter_class=argparse.RawTextHelpFormatter,
                                        usage="%s  [-h] [-p pluginName [binVersion ...]]" %
                                              (' '.join(args[:2])),
                                        epilog="Example: %s -p ctffind -p eman 2.11  \n\n%s" %
                                               (' '.join(args[:2]), pluginRepo.printPluginInfo()))
uninstallParser.add_argument('-p', '--plugin', action='append', nargs='+',
                             metavar=('pluginName', 'binVersion'),
                             help='- pluginName: Name of the plugin to uninstall from the list of\n'
                                  '              installed plugins shown below. \n'
                                  '- binVersion: optional param to control the binary version uninstalled\n'
                                  "              - If none given, will uninstall all binaries installed and\n"
                                  "                the plugin's python package. (See ctffind example below)\n"
                                  "              - If binVersion specified, it will only uninstall this one, \n"
                                  "                leaving the python module intact.(See eman example below)\n"
                                  'Several plugins may be specified in a single command, \n'
                                  'each following the "-p pluginName [binVersion]" pattern. \n')


parsedArgs = parser.parse_args(args[1:])
# if parsedArgs.mode == MODE_INSTALL_PLUGIN:
for cmdTarget in parsedArgs.plugin:
    pluginName = cmdTarget[0]
    plugin = pluginRepo.pluginDict.get(pluginName, None)
    if plugin:
        binVersion = "" if len(cmdTarget) <= 1 else cmdTarget[1]
        if parsedArgs.mode == MODE_INSTALL_PLUGIN:
            plugin.installPipModule()
            if not parsedArgs.noBin and binVersion != '0':
                binArgs = [Environment._getExtName(pluginName, binVersion)] if binVersion else []
                plugin.installBin(args=binArgs)
        elif parsedArgs.mode == MODE_UNINSTALL_PLUGIN:
            plugin.uninstallBin(version=binVersion)
            if not binVersion:  # only remove pip plugin if no version specified
                plugin.uninstallPip()
    else:
        print("Can't (un)install %s plugin because it is not installed or doesn't exist." % pluginName)

# elif parsedArgs.mode == MODE_UNINSTALL_PLUGIN:
#     for cmdTarget in parsedArgs.plugin:
#         pluginName = cmdTarget[0]
#         plugin = pluginRepo.pluginDict.get(pluginName, None)
#         if plugin:
#             binVersion = "" if len(cmdTarget) <= 1 else cmdTarget[1]
#             plugin.uninstallBin(version=binVersion)
#             if not binVersion: # only remove pip plugin if no version specified
#                 plugin.uninstallPip()
#         else:
#             print("Can't uninstall %s plugin because it is not installed or doesn't exist.")



