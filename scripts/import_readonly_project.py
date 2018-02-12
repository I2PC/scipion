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

import sys, os
from pyworkflow.manager import Manager
from pyworkflow.project import PROJECT_RUNS, PROJECT_DBNAME
import pyworkflow.utils as pwutils


def usage(error):
    print """
    ERROR: %s

    Usage: scipion python scripts/import_readonly_project.py FROM_PATH [PROJECT_NAME] [LOCATION]
    FROM_PATH: Path to the read only project we want to import
    PROJECT_NAME (optional): A name for our project. If none given, will use the one in the read only project
    LOCATION (optional): where to create our hybrid project.  Defaults to scipion default project location.
    e.g: scipion python scripts/create_readonly_copy.py  /path/to/readonly/project MyHybridProject /some/path

    """ % error
    sys.exit(1)

n = len(sys.argv)

if n < 2 or n > 5:
    usage("Incorrect number of input parameters")
fromLocation = sys.argv[1]
projName = fromLocation.rsplit('/', 1)[-1] if n < 3 else sys.argv[2]

jsonFile = None if n < 4 else os.path.abspath(sys.argv[3])
location = None if n < 5 else sys.argv[4]

path = os.path.join(os.environ['SCIPION_HOME'], 'pyworkflow', 'gui', 'no-tkinter')
sys.path.insert(1, path)

# Create a new project
manager = Manager()

if manager.hasProject(projName):
    usage("There is already a project with this name: %s"
          % pwutils.red(projName))

project = manager.createProject(projName, location=location)
# copy sqlite and make links of runs
pwutils.copyFile(os.path.join(fromLocation, PROJECT_DBNAME),
                 project.path)
fromRuns = os.path.join(fromLocation, PROJECT_RUNS)
for r in os.listdir(fromRuns):
    pwutils.createAbsLink(os.path.join(fromRuns, r),
                          os.path.join(project.runsPath, r))