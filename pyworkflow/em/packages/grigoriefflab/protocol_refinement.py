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
from pyworkflow.utils import join
import pyworkflow.em as em 
from protocol_frealign_base import ProtFrealignBase
from convert import readSetOfParticles


class ProtFrealign(ProtFrealignBase, em.ProtRefine3D):
    """ Protocol to refine a 3D map using Frealign. The algorithms implemented
are optimized to perform  efficiently the correction for the contrast
transfer function of the microscope and refinement of three-dimensional
reconstructions.
    """
    _label = 'frealign'
    
    def __init__(self, **args):
        ProtFrealignBase.__init__(self, **args)
    
    def createOutputStep(self):
        lastIter = self._getLastIter()
        inputSet = self._getInputParticles()

        # Register output volume
        volFn = self._getFileName('iter_vol', iter=lastIter)
        vol = em.Volume()
        vol.setFileName(volFn)
        vol.setSamplingRate(inputSet.getSamplingRate())
        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputParticles, vol)

        # Register output Particles with their 3D alignment
        file2 = self._getFileName('output_par', iter=lastIter)
        partSet = self._createSetOfParticles()
        partSet.copyInfo(inputSet)
        readSetOfParticles(inputSet, partSet, file2)
        self._defineOutputs(outputParticles=partSet)
        self._defineTransformRelation(self._getInputParticlesPointer(), partSet)
        
        if not self.doContinue:
            self._defineSourceRelation(self.input3DReference, vol)
            self._defineSourceRelation(self.input3DReference, partSet)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _citations(self):
        return ['Lyumkis2013', 'Sindelar2012', 'Grigorieff2007', 'Wolf2006', 'Stewart2004', 'Grigorieff1998']
    
