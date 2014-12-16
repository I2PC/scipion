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
Update project settings using the current settings from
ProjectSettings().loadConfig().
"""

import sys
from os.path import join

from pyworkflow.utils.path import copyFile, cleanPath
from pyworkflow.config import ProjectSettings, pw
from pyworkflow.manager import Manager


def updateSettings():
    """ Write the settings.sqlite configuration file for all projects.
    """
    # Load the default settings and write them in pw.SETTINGS (a sqlite file).
    settings = ProjectSettings()
    settings.loadConfig()
    settings.write(pw.SETTINGS)

    # Update the settings in all existing projects, by copying the sqlite file.
    manager = Manager()
    for p in manager.listProjects():
        fpath = join(manager.getProjectPath(p.getName()), 'settings.sqlite')
        print "Copying settings to:", fpath
        copyFile(pw.SETTINGS, fpath)

    cleanPath(pw.SETTINGS) # we dont need to store the settings.sqlite file



if __name__ == '__main__':
    if len(sys.argv) > 1:
        sys.exit("Error: %s does not take any argument (got: %s)." %
                 sys.argv[0], " ".join(sys.argv[1:]))
    updateSettings()
