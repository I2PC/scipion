#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     Yaiza Rancel (yrancel@cnb.csic.es)
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

import sys
import os
from time import sleep

import pyworkflow as pw
from pyworkflow.project import Manager
from pyworkflow.project import Project
import pyworkflow.utils as pwutils


def usage(error):
    print """
    ERROR: %s

    Usage: scipion python scripts/refresh_project.py project_name [refresh_period]
        This script will refresh the project.sqlite file of the specified project.
        We can specify how often it will refresh in seconds (defaults 60).
        e.g.
        scipion python scripts/refresh_project.py TestExtractParticles 15
    """ % error
    sys.exit(1)


n = len(sys.argv)

if n > 3:
    usage("This script accepts 2 parameters: the project name (mandatory) and "
          "the refresh period in seconds (optional).")

projName = sys.argv[1]
if n == 3:
    wait = int(sys.argv[2])
else:
    wait = 60

path = pw.join('gui', 'no-tkinter')
sys.path.insert(1, path)

# Create a new project
manager = Manager()

if not manager.hasProject(projName):
    usage("There is no project with this name: %s"
          % pwutils.red(projName))

# the project may be a soft link which may be unavailable to the cluster so get the real path
try:
    projectPath = os.readlink(manager.getProjectPath(projName))
except:
    projectPath = manager.getProjectPath(projName)

project = Project(projectPath)
project.load()
while True:
    runs = project.getRuns()
    sleep(wait)

