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
This sub-package will contains Spider protocols
"""

from spider import *

from protocol_filters import SpiderProtFilter
from protocol_align_apsr import SpiderProtAlignAPSR
from protocol_custommask import SpiderProtCustomMask
from protocol_ca_pca import SpiderProtCAPCA
from protocol_ward import SpiderProtClassifyWard
from protocol_mda import SpiderWfMDA
from protocol_align_pairwise import SpiderProtAlignPairwise

from wizard import SpiderProtMaskWizard, SpiderProtMaskRadiiWizard, SpiderGaussianFilterWizard #, SpiderFermiFilterWizard

from viewer import SpiderViewerGeneric
from viewer_capca import SpiderViewerCAPCA
from viewer_ward import SpiderViewerWard
