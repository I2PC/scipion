# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package contains data and protocol classes
wrapping Localized recontruction of subunits.
"""
from bibtex import _bibtex # Load bibtex dict with references

_logo = "opic_logo.png"
LOCALREC_HOME = 'LOCALREC_HOME'
RELION_HOME = 'RELION_HOME'
TMP_RELION_HOME = 'TMP_RELION_HOME'
LOCALREC_RELION_HOME = 'LOCALREC_RELION_HOME'

from convert import *
from protocol_localized import ProtLocalizedRecons
from protocol_localized_extraction import ProtLocalizedExtraction


def validateInstallation():
    """ This function will be used to check if package is properly installed."""
    missingPaths = ["%s: %s" % (var, os.environ[var])
                    for var in [LOCALREC_HOME, LOCALREC_RELION_HOME]
                    if not os.path.exists(os.environ[var])]

    if missingPaths:
        return ["Missing variables:"] + missingPaths
    else:
        return [] # No errors
