#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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


def usage(error):
    print """
    ERROR: %s
    
    Usage: scipion python scripts/clean_projects.py [SCIPION_USER_DATA] [--delete]
        Clean projects that has expired (now - creation time > life time)
        Optional to pass SCIPION_USER_DATA folder from which to read 'projects'.
        If --delete is not passed, only a message with the projects to be deleted
        will be shown. If used --delete, the projects will be deleted from filesystem.
    """ % error
    sys.exit(1)    

n = len(sys.argv)

if n > 3:
    usage("Incorrect number of input parameters")
    
delete = '--delete' in sys.argv

arg1 = sys.argv[1]
if n > 1 and arg1 != '--delete':
    customUserData = arg1
else:
    customUserData = os.environ['SCIPION_USER_DATA']

print "Loading projects from:\n", customUserData 
 
# Create a new project
manager = Manager(SCIPION_USER_DATA=customUserData)

for projInfo in manager.listProjects():
    projName = projInfo.getName()
    proj = manager.loadProject(projName)
    settings = proj.getSettings()
    
    leftTime = proj.getLeftTime()
    if (leftTime is not None and 
        leftTime.days < 0):
        if delete:
            print "Deleting: %s (%s) " % (projName, leftTime)
            manager.deleteProject(projName)
        else: 
            print "Should delete: %s (%s)" % (projName, leftTime)
    
    

