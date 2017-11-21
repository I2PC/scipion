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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package contains Relion protocols and tools.
"""

from bibtex import _bibtex # Load bibtex dict with references
_logo = "relion_logo.png"
_references = ['Scheres2012a', 'Scheres2012b', 'Chen2012']

from convert import *
from protocol_create_mask3d import ProtRelionCreateMask3D
from protocol_classify2d import ProtRelionClassify2D
from protocol_classify3d import ProtRelionClassify3D
from protocol_refine3d import ProtRelionRefine3D
from protocol_reconstruct import ProtRelionReconstruct
from protocol_postprocess import ProtRelionPostprocess
from protocol_preprocess import ProtRelionPreprocessParticles
from protocol_polish import ProtRelionPolish
from protocol_sort import ProtRelionSortParticles
from protocol_subtract import ProtRelionSubtract
from protocol_expand_symmetry import ProtRelionExpandSymmetry
from protocol_initialmodel import ProtRelionInitialModel
from protocol_localres import ProtRelionLocalRes

from protocol_autopick import ProtRelionAutopickFom, ProtRelionAutopick
from protocol_autopick_v2 import ProtRelion2Autopick
from protocol_extract_particles import ProtRelionExtractParticles

from protocol_export_ctf import ProtRelionExportCtf

# Wizards
from wizard import *

from viewer import *

_environ = getEnviron()


def validateInstallation():
    """ This function will be used to check if RELION is properly installed. """
    missingPaths = ["%s: %s" % (var, _environ[var])
                    for var in [RELION_HOME]
                    if not os.path.exists(_environ[var])]

    if missingPaths:
        return ["Missing variables:"] + missingPaths
    else:
        return [] # No errors
