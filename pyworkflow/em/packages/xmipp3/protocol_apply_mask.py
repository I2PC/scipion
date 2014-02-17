# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'xmipp@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package contains protocols for applying 3D masks.
"""

from pyworkflow.em import *
from geometrical_mask import XmippGeometricalMask3D

SOURCE_GEOMETRY=0
SOURCE_MASK=1

class XmippProtApplyMask3D(ProtPreprocessVolumes,XmippGeometricalMask3D):
    """ Apply a 3D mask to a volume or set of volumes """
    _label = 'apply mask volumes'
    
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('volume', PointerParam, pointerClass="SetOfVolumes", label="Input volume",
                      help="Volume or set of volumes to be masked")
        form.addParam('source', EnumParam, label='Mask source', default=SOURCE_GEOMETRY, choices=['Geometry','Created mask'])
        XmippGeometricalMask3D.defineParams(self, form, isGeometry='source==%d'%SOURCE_GEOMETRY, addSize=False)
        form.addParam('inputMask', PointerParam, pointerClass="VolumeMask", label="Input mask",condition='source==%d'%SOURCE_MASK)

    def _insertAllSteps(self):
        self.maskedVolumes = self._getPath('volumes.stk')
        self._insertFunctionStep('copyVolumes')
        self._insertFunctionStep('applyMask')
        self._insertFunctionStep('createOutput')
    
    def copyVolumes(self):
        # COSS Hace falta manejar el caso de que sea un volumen o varios
        self.volume.get().writeStack(self.maskedVolumes)

    def applyMask(self):
        args="-i %s "%self.maskedVolumes
        if self.source.get()==SOURCE_GEOMETRY:
            (Xdim, _, _, self.ndim)=self.volume.get().getDimensions()
            args += XmippGeometricalMask3D.argsForTransformMask(self,Xdim)
        elif self.source.get()==SOURCE_MASK:
            args+="--mask binary_file %s"%self.inputMask.get().getFileName()
        self.runJob("xmipp_transform_mask",args)

    def createOutput(self):
        maskedVolumes = self._createSetOfVolumes()
        maskedVolumes.setSamplingRate(self.volume.get().getSamplingRate())
        for i in range(1,self.ndim+1):
            vol=Volume()
            vol.setLocation(i,self.maskedVolumes)
            maskedVolumes.append(vol)
        self._defineOutputs(maskedVolumes=maskedVolumes)
        self._defineTransformRelation(self.volume.get(), maskedVolumes)
        
    def _summary(self):
        messages = []      
        if self.source.get()==SOURCE_GEOMETRY:
            messages+=XmippGeometricalMask3D.summary(self)
        return messages
    