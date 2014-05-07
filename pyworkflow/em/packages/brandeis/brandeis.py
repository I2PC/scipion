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
This sub-package contains data and protocol classes
wrapping Grigrorieff Lab programs at Brandeis
"""

from os.path import join
from pyworkflow.em import *
from constants import *

CTFFIND = 'ctffind3.exe'
FREALIGN = 'frealign_v9.exe'
FREALIGNMP = 'frealign_v9_mp.exe'
CALC_OCC = 'calc_occ.exe'
RSAMPLE = 'rsample.exe'

CTFFIND_PATH = join(os.environ['CTFFIND_HOME'], CTFFIND)
FREALIGN_PATH = join(os.environ['FREALIGN_HOME'], 'bin', FREALIGN)
FREALIGNMP_PATH = join(os.environ['FREALIGN_HOME'], 'bin', FREALIGNMP)
CALC_OCC_PATH = join(os.environ['FREALIGN_HOME'], 'bin', CALC_OCC)
RSAMPLE_PATH = join(os.environ['FREALIGN_HOME'], 'bin', RSAMPLE)