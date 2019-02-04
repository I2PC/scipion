#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
"""
Launch main project window 
"""

import sys
import os

import time

from pyworkflow.project import Manager
from pyworkflow.gui.project import ProjectWindow
import pyworkflow.utils as pwutils


if __name__ == '__main__':

    # Add callback for remote debugging if available.
    try:
        from rpdb2 import start_embedded_debugger
        from signal import signal, SIGUSR2
        signal(SIGUSR2, lambda sig, frame: start_embedded_debugger('a'))
    except ImportError:
        pass

    if len(sys.argv) > 1:
        manager = Manager()
        projName = os.path.basename(sys.argv[1])

        # Handle special name 'here' to create a project
        # from the current directory
        if projName == 'here':
            cwd = os.environ['SCIPION_CWD']
            print "\nYou are trying to create a project here:", pwutils.cyan(cwd)

            if os.listdir(cwd):
                print pwutils.red('\nWARNING: this folder is not empty!!!')
            key = raw_input("\nDo you want to create a project here? [y/N]?")

            if key.lower().strip() != 'y':
                print "\nAborting..."
                sys.exit(0)
            else:
                print "\nCreating project...."
                projName = os.path.basename(cwd)
                projDir = os.path.dirname(cwd)
                proj = manager.createProject(projName, location=projDir)

        elif projName == 'last':  # Get last project
            projects = manager.listProjects()
            if not projects:
                sys.exit("No projects yet, cannot open the last one.")
            projName = projects[0].projName
            
        projPath = manager.getProjectPath(projName)
        # try:
        projWindow = ProjectWindow(projPath)
        # except Exception as e:
        #     # Print any exception
        #     print("ERROR: At pw_project.py loading Project %s.\n"
        #           "       Message: %s\n" % (projPath, e))
        #
        #     import traceback
        #     traceback.print_exc(file=sys.stderr)
        #
        #     sys.exit(e)

        projWindow.show()
    else:
        print "usage: pw_project.py PROJECT_NAME"
