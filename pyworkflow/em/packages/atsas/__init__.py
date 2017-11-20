# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
This sub-package contains data and protocol classes
wrapping ATSAS programs http://www.embl-hamburg.de/biosaxs/software.html
"""
from pyworkflow.utils import commandExists

_logo = "atsas_logo.gif"
CRYSOL = "crysol"
from bibtex import _bibtex # Load bibtex dict with references
from atsas import *
from protocol_pdb_to_saxs import AtsasProtConvertPdbToSAXS
from viewer import AtsasViewer

def validateInstallation():
    """ This function will be used to check if ATSAS is properly installed. """
    missingPaths = []

    if not (commandExists(CRYSOL)):
        missingPaths.append("%s not found in the system" % CRYSOL)

    if missingPaths:
        return ["Missing variables:"] + missingPaths
    else:
        return [] # No errors
