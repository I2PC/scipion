# **************************************************************************
# *
# * Authors:     Airen Zaldivar Peraza (azaldivar@cnb.csic.es)
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
In this module are protocol base classes related to EM.
There should be sub-classes in the different packages from
each EM-software package.
"""
from protocol import *
from protocol_import import *
from protocol_micrographs import *
from protocol_movies import *
from protocol_particles import *
from protocol_2d import *
from protocol_3d import *
from protocol_sets import *
from protocol_tiltpairs import *
from protocol_ctf_assign import ProtCTFAssign
from protocol_alignment_assign import ProtAlignmentAssign
from protocol_batch import *
from protocol_classes_consensus import ProtClassesConsensus, ViewerClassesConsensus

from parallel import ProtTestParallel

