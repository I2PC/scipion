# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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
This EM module contains Gautomatch auto-picking protocol 
"""

_logo = "gautomatch_logo.png"
GAUTOMATCH_HOME = 'GAUTOMATCH_HOME'

from convert import getEnviron
from bibtex import _bibtex # Load bibtex dict with references

from protocol_gautomatch import ProtGautomatch
from viewer import GautomatchViewer
from wizard import *
_environ = getEnviron()

def validateInstallation():
    """ This function will be used to check if package is properly installed."""
    missingPaths = ["%s: %s" % (var, _environ[var])
                    for var in [GAUTOMATCH_HOME]
                    if not os.path.exists(_environ[var])]

    if missingPaths:
        return ["Missing variables:"] + missingPaths
    else:
        return [] # No errors
