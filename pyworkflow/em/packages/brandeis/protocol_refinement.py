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
    """ This class implements the wrapper to single particle refinement protocol with frealign."""
    _label = 'frealign'
    _references = ['[[http://dx.doi.org/10.1016/j.jsb.2006.05.004][Grigorieff N,  JSB (2007)]]',
                   '[[http://www.ncbi.nlm.nih.gov/pubmed/16384646][Wolf M, et.al, Ultramicroscopy (2006)]]',
                   '[[http://www.ncbi.nlm.nih.gov/pubmed/15556702][Stewart A & Grigorieff N, Ultramicroscopy (2004)]]',
                   '[[http://www.ncbi.nlm.nih.gov/pubmed/9571020][Grigorieff N, JMB (1998)]]',
                   '[[http://www.sciencedirect.com/science/article/pii/S104784771200144X][Sindelar CV & Grigorieff N, Ultramicroscopy (2012)]]',
                   '[[http://www.sciencedirect.com/science/article/pii/S1047847713001858][Lyumkis D, et. al, JSB (2013)]]'
                    ]
    
    def __init__(self, **args):
        ProtFrealignBase.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL