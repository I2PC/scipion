#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
from datetime import datetime

from pyworkflow.manager import Manager
from pyworkflow.gui.project import ProjectWindow
import pyworkflow.utils as pwutils


def usage(error):
    print """
    ERROR: %s
    
    Usage: scipion python scripts/edit_workflow.py workflow.json
        Edit the provide json file with scipion workflow.
        It will create a project, import the workflow and save
        the workflow back before closing the project.
        After that, the project will be deleted.
    """ % error
    sys.exit(1)    

n = len(sys.argv)

if n < 2 or n > 3:
    usage("Incorrect number of input parameters")
    
jsonFn = os.path.abspath(sys.argv[1])

now = datetime.now()
tempSpace = "editor-%s" % now.strftime('%Y%m%d-%H%M%S')
customUserData = os.path.join(os.environ['SCIPION_USER_DATA'],
                              'tmp', tempSpace)

pwutils.makePath(os.path.join(customUserData, 'projects'))

print "Loading projects from:\n", customUserData 
 
# Create a new project
manager = Manager(SCIPION_USER_DATA=customUserData)

projName = os.path.basename(jsonFn)
proj = manager.createProject(projName)
projPath = manager.getProjectPath(projName)
proj.loadProtocols(jsonFn)

class EditorProjectWindow(ProjectWindow):
    def close(self, e=None):
        try:
            print "Writing protocols to: ", jsonFn
            proj.getRunsGraph() # Build project runs graph
            proj.exportProtocols(proj.getRuns(), jsonFn)
            print "Deleting temporary folder: ", customUserData
            pwutils.cleanPath(customUserData)
        except Exception, ex:
            print "Error saving the workflow: ", ex
        ProjectWindow.close(self, e)

projWindow = EditorProjectWindow(projPath)
projWindow.show()


