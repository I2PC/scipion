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
This modules serve to define some Configuration classes
mainly for project GUI
"""

from pyworkflow.utils.path import cleanPath
from pyworkflow.config import *


def writeDefaults():
    settings = ProjectSettings()
    settings.loadConfig()

    #writeConfig(config, 'configuration.xml')
    dbPath = pw.SETTINGS
    print "Writing default settings to: ", dbPath
    if exists(dbPath):
        dbPathBackup = '%s.%d' % (dbPath, int(time.time()))
        print "File %s exists, moving to %s ..." % (dbPath, dbPathBackup)
        os.rename(dbPath, dbPathBackup)
    settings.write(dbPath)


def updateSettings():
    """ Write the settings.sqlite default configuration
    and also update each project settings.
    """
    writeDefaults()
    # Update the settings to all existing projects
    from pyworkflow.manager import Manager
    from pyworkflow.utils.path import copyFile

    manager = Manager()
    projects = manager.listProjects()

    for p in projects:
        proj = manager.loadProject(p.getName())
        projSettings = proj.settingsPath
        print "Copying settings to: ", join(p.getName(), projSettings)
        copyFile(pw.SETTINGS, projSettings)
        
    cleanPath(pw.SETTINGS) # we dont need to store the settings.sqlite 



if __name__ == '__main__':
    #Write default configurations
    updateSettings()
