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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains the protocol to obtain a refined 3D recontruction from a set of particles using Frealign
"""
import os
from pyworkflow.utils import *
from pyworkflow.em import *
from data import *
from brandeis import *
from constants import *
from protocol_frealign_base import ProtFrealignBase


class ProtFrealign(ProtFrealignBase, ProtRefine3D):
    """ Protocol to refine a 3D map using Frealign. The algorithms implemented
are optimized to perform  efficiently the correction for the contrast
transfer function of the microscope and refinement of three-dimensional
reconstructions.
    """
    _label = 'frealign'

    
    def __init__(self, **args):
        ProtFrealignBase.__init__(self, **args)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _citations(self):
        return ['Lyumkis2013', 'Sindelar2012', 'Grigorieff2007', 'Wolf2006', 'Stewart2004', 'Grigorieff1998']
    
