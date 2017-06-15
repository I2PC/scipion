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
from pyworkflow.object import Boolean


def usage(error):
    print """
    ERROR: %s
    
    Usage: scipion python scripts/config_projects.py ProjectName [readOnly=True|False] [lifeTime=X|None]
        This script show (or edit) some of the configuration of the project.
        Use readOnly=True (or False) to set/unset read only property
        Use lifeTime=X for setting X hours or None to unset life time of the project.
    """ % error
    sys.exit(1)    

n = len(sys.argv)

if n < 2 or n > 4:
    usage("Incorrect number of input parameters")
    
# Load the given project
projectsDir = os.path.join(os.environ['SCIPION_USER_DATA'], 'projects')
projName = sys.argv[1]
manager = Manager()
project = manager.loadProject(projName)

if project is None:
    usage("Project '%s' does not exist in: \n  %s" % (projName, projectsDir))

setReadOnly = False
setLifeTime = False

for arg in sys.argv:
    if arg.startswith('readOnly='):
        setReadOnly = True
        value = arg.split('readOnly=')[1]
        b = Boolean(value=value)
        readOnlyValue = b.get()
    elif arg.startswith('lifeTime='):
        setLifeTime = True
        value = arg.split('lifeTime=')[1]
        lifeTimeValue = None if value == 'None' else int(value)

if setReadOnly:
    project.setReadOnly(readOnlyValue)
    
if setLifeTime:
    project.settings.setLifeTime(lifeTimeValue)
    
if setReadOnly or setLifeTime:
    # Truly write settings
    project.settings.write()

print "Projects: ", projectsDir
print "Project name: ", projName
print " Settings: "
print "   readOnly = %s" % project.isReadOnly()
print "   lifeTime = %s (hours)" % project.settings.getLifeTime()
   

