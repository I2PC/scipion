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
This sub-package contains data and protocol classes
wrapping Grigorieff Lab programs at Brandeis
"""
from bibtex import _bibtex # Load bibtex dict with references

_logo = "brandeis_logo.png"

from grigoriefflab import *

from protocol_ctffind import ProtCTFFind
from protocol_refinement import ProtFrealign
from protocol_magdist_estimate import ProtMagDistEst
from protocol_magdist_correct import ProtMagDistCorr
from protocol_magdist_correct_coords import ProtMagDistCorrCoord
from protocol_ml_classification import ProtFrealignClassify
from protocol_unblur import ProtUnblur
from protocol_summovie import ProtSummovie
from viewer import FrealignViewer, MagDistEstViewer
# Wizards
from wizard import *
