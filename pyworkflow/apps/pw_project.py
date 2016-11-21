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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
Launch main project window 
"""

import sys
import os

from pyworkflow.manager import Manager
from pyworkflow.gui.project import ProjectWindow


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
        if projName == 'last':  # Get last project
            projects = manager.listProjects()
            if not projects:
                sys.exit("No projects yet, cannot open the last one.")
            projName = projects[0].projName
            
        projPath = manager.getProjectPath(projName)
        try:
            projWindow = ProjectWindow(projPath)
        except Exception as e:
            sys.exit(e)
        projWindow.show()
    else:
        print "usage: pw_project.py PROJECT_NAME"
