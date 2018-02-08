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
This package contains the protocols and data for EMAN2
"""

from bibtex import _bibtex # Load bibtex dict with references
from pyworkflow.plugin import Plugin
import os

EMAN_DIR_VAR = 'EMAN2DIR'
VARS = {EMAN_DIR_VAR: '/sample/plugin/default/path/software/em/eman-2.12'}

_plugin = Plugin('eman2',
                 version=2,
                 configVars=VARS,
                 logo="eman2_logo2.png",
                 references=['Tang2007'])

from eman2 import *
from protocol_boxing import EmanProtBoxing
from protocol_initialmodel import EmanProtInitModel
from protocol_reconstruct import EmanProtReconstruct
from protocol_refineasy import EmanProtRefine
from protocol_autopick import SparxGaussianProtPicking
from viewer import EmanViewer, RefineEasyViewer
from wizard import SparxGaussianPickerWizard


_logo = _plugin.logo
_references = _plugin.references
_environ = getEnviron()

def validateInstallation():
    """ This function will be used to check if package is properly installed."""
    missingPaths = ["%s: %s" % (var, _environ[var])
                    for var in VARS
                    if not os.path.exists(_environ[var])]

    if missingPaths:
        return ["Missing variables:"] + missingPaths
    else:
        return [] # No errors
