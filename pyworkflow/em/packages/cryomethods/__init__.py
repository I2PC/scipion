# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *              Javier Vargas Balbuena (javier.vargasbalbuena@mcgill.ca)
# *
# * Department of Anatomy and Cell Biology, McGill University
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
This sub-package contains cryoMethods protocols and tools.
"""

# from bibtex import _bibtex # Load bibtex dict with references
_logo = "cryomethods_logo.png"
_references = []
import os
from convert import getSupportedVersions, getVersion, getEnviron
from initial_volume_selector import ProtInitialVolumeSelector
#
# # Wizards
# from wizard import *
from viewer import *
#
_environ = getEnviron()


def validateInstallation():
    """ This function will be used to check if RELION is properly installed. """
    missingPaths = ["%s: %s" % (var, _environ[var])
                    for var in ['RELION_HOME']
                    if not os.path.exists(_environ[var])]

    if missingPaths:
        return ["Missing variables:"] + missingPaths
    else:
        return [] # No errors
