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
This sub-package contains data and protocol classes
wrapping Kai Zhang's GCTF program
"""
import os

_logo = "gctf_logo.png"
GCTF_HOME = 'GCTF_HOME'

from bibtex import _bibtex # Load bibtex dict with references
from convert import getEnviron

from protocol_gctf import ProtGctf
from protocol_gctf_refine import ProtGctfRefine

# Wizards
from wizard import GctfCTFWizard
_environ = getEnviron()

# We need this import to register the specific viewing command
# when visualizing Gctf results
from viewer import GctfViewer


def validateInstallation():
    """ This function will be used to check if package is properly installed."""
    missingPaths = ["%s: %s" % (var, _environ[var])
                    for var in [GCTF_HOME]
                    if not os.path.exists(_environ[var])]

    if missingPaths:
        return ["Missing variables:"] + missingPaths
    else:
        return [] # No errors
