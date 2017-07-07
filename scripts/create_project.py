#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
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

import sys, os
from pyworkflow.manager import Manager
import pyworkflow.utils as pwutils


def usage(error):
    print """
    ERROR: %s

    Usage: scipion python scripts/create_project.py
        name="project name"
        [workflow="file"] path to a Scipion json workflow
        [location="folder"] where to create it, defaults to scipion default location
        This script will create a project project, optionally based on a workflow file
    """ % error
    sys.exit(1)


n = len(sys.argv)

if n < 2 or n > 4:
    usage("Incorrect number of input parameters")

projName = sys.argv[1]

jsonFile = None if n < 3 else os.path.abspath(sys.argv[2])
location = None if n < 4 else sys.argv[3]

path = os.path.join(os.environ['SCIPION_HOME'], 'pyworkflow', 'gui', 'no-tkinter')
sys.path.insert(1, path)

# Create a new project
manager = Manager()

if manager.hasProject(projName):
    usage("There is already a project with this name: %s"
          % pwutils.red(projName))

if jsonFile is not None and not os.path.exists(jsonFile):
    usage("Unexistent json file: %s" % pwutils.red(jsonFile))

project = manager.createProject(projName, location=location)

if jsonFile is not None:
    protDict = project.loadProtocols(jsonFile)