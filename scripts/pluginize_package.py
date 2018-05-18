#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     Yaiza Rancel (cyrancel@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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
import sys
import pyworkflow.utils as pwutils
from pyworkflow.em import PACKAGES_PATH, getPackages
from pyworkflow import findResource
import argparse

pipEmptyFiles = ['CHANGES.txt', 'MANIFEST.in', 'README.rst']
templatesFolder = os.path.join(os.environ['SCIPION_HOME'], 'config', 'templates')
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 usage='Script to create basic plugin structure for an existing Scipion package \n'
                                       'in an empty folder (which will become the plugin repository)\n'
                                       'and make a link to it in the original location (pyworkflow/em/packages).\n',
                                 epilog='Examples: \n'
                                        'scipion python scripts/pluginize_package.py grigoriefflab '
                                        '~/scipion_grigoriefflab\n'
                                        'scipion python scripts/pluginize_package.py --help\n')
parser.add_argument('package_name', help='Name of the package to pluginize')
parser.add_argument('plugin_dir', help='Empty directory outside of Scipion where we will move the plugin. '
                                       "Will be created if it doesn't exist")
if len(sys.argv) < 2:
    print("ERROR: Wrong number of arguments!!!")
    parser.print_help()
    sys.exit(0)

args = parser.parse_args()

pluginDir = args.plugin_dir
packageName = args.package_name
allPackages = getPackages()
package = allPackages.get(packageName, None)

if package is None:
    print("WARNING: Package %s does not exist" % packageName)
    print('The list of available packages is:\n%s' % '\n'.join(allPackages.keys()))
    parser.print_help()
else:
    # 1. Make new dir
    if not pwutils.exists(pluginDir):
        print('Creating new folder for plugin in %s' % pluginDir)
        pwutils.makePath(pluginDir)
    elif len(os.listdir(pluginDir)) != 0:
        print("Directory %s is not empty. Please enter an empty or non existing directory" % pluginDir)
        sys.exit(0)

    # 2. Move logo to package folder
    logo = getattr(package, '_logo', None)
    if logo:
        logoPath = findResource(logo)
        print('Found logo %s. Moving it to plugin folder....' % logo)
        pwutils.moveFile(logoPath, os.path.join(PACKAGES_PATH, packageName, logo))

    # 3. Move everything to new dir created in 1
    pwutils.moveTree(os.path.join(PACKAGES_PATH, packageName), os.path.join(pluginDir, packageName))

    # 4. Make link to new dir
    pwutils.createAbsLink(os.path.join(pluginDir, packageName),
                          os.path.join(PACKAGES_PATH, packageName))

    # 5. Copy plugin files
    pwutils.copyFile(os.path.join(templatesFolder, 'setup.template'),
                     os.path.join(pluginDir, 'setup.py'))
    pwutils.copyFile(os.path.join(templatesFolder, 'plugin.template'),
                     os.path.join(pluginDir, packageName, 'plugin.py'))

    # 6. Copy pip files
    for f in pipEmptyFiles:
        open(os.path.join(pluginDir, f), 'a').close()


    # 7. Check if requests is installed
    try:
        import requests
    except ImportError as e:
        print('Requests pip package not found. Installing...')
        os.system('%s install requests --no-xmipp' % os.path.join(os.environ['SCIPION_HOME'], 'scipion'))

