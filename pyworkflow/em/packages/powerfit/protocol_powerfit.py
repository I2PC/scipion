# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), May 2013
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

import math
from glob import glob
import numpy

from pyworkflow.em import *  
from pyworkflow.utils import * 
from pyworkflow.em.convert import ImageHandler
from pyworkflow.protocol.constants import LEVEL_ADVANCED
import powerfit
from pyworkflow.em.packages.ccp4.convert import Ccp4Header

# TODO: Move to 3D Tools
class PowerfitProtRigidFit(ProtFitting3D):
    """ Protocol for fitting a PDB into a 3D volume
    
    This is actually a wrapper to the program Powerfit.
    See documentation at:
       http://www.bonvinlab.org/education/powerfit
    """
    _label = 'rigid fit'
    
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVol', PointerParam, pointerClass='Volume',
                      label="Input volume", important=True)
        form.addParam('inputPDB', PointerParam, pointerClass='PdbFile',
                      label="Input PDB", important=True)
        form.addParam('resolution', FloatParam, default=6,
                      label="Resolution (A)", important=True, help="Resolution for the fitting. The PDB is filtered to this frequency.")
        form.addParam('angleStep',FloatParam, label="Angular step", default=10.0, expertLevel=LEVEL_ADVANCED,
                      help='Angular step for the alignment search')
        form.addParam('doLaplacian',BooleanParam, label="Apply Laplacian", default=False, expertLevel=LEVEL_ADVANCED,
                      help='Apply a Laplacian to the volume to highlight borders')
        form.addParam('doCoreWeight',BooleanParam, label="Apply core weight", default=False, expertLevel=LEVEL_ADVANCED,
                      help='Apply core weights')
        form.addParam('otherPowerfit', StringParam, default='', expertLevel=LEVEL_ADVANCED,
                      label='Other parameters for Powerfit',
                      help='See http://www.bonvinlab.org/education/powerfit') 
        form.addParallelSection(threads=4, mpi=0)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('powerfitWrapper')
        
    #--------------------------- STEPS functions -------------------------------
    def powerfitWrapper(self):
        img = ImageHandler()
        fnVol = self._getExtraPath('volume.mrc')
        vol = self.inputVol.get()
        img.convert(vol,fnVol)

        ccp4header = Ccp4Header(fnVol, readHeader= True)
        ccp4header.setOffset(vol.getOrigin())
        ccp4header.setSampling(vol.getSamplingRate())
        ccp4header.writeHeader()

        args = "%s %f %s -d %s -p %d -a %f"%(fnVol,self.resolution,self.inputPDB.get().getFileName(), self._getExtraPath(),self.numberOfThreads,
                                             self.angleStep)
        if self.doLaplacian:
            args+=" -l"
        if self.doCoreWeight:
            args+=" -cw"
        self.runJob("powerfit",args)
        
    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        summary.append('Angular step: %d' % self.angleStep)
        if self.doLaplacian:
            summary.append("Apply Laplacian")
        if self.doLaplacian:
            summary.append("Apply core weights")
        if not self.otherPowerfit.empty():
            summary.append('Other powerfit parameters: %s' % self.otherPowerfit)
        return summary

    def _methods(self):
        summary = []
        summary.append('We rigidly fitted the structure %s into the volume with an angular step of %f using Powerfit [vanZundert2015].'
                       % (self.inputPDB.get().getNameId(),self.inputVol.get().getNameId(),self.angleStep))
        if self.doLaplacian:
            summary.append("We applied a Laplacian filter to the input volume.")
        if self.doLaplacian:
            summary.append("We used core weighted local cross-correlations.")
        return summary

    def _citations(self):
        return ['vanZundert2015']
    
    def _validate(self):
        errors = []
        if which('powerfit') is '':
            errors.append('You should have the program powerfit in the PATH')
        return errors
        