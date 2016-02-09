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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.protocol.params import MultiPointerParam, PointerParam, BooleanParam
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.packages.xmipp3.convert import getImageLocation
import os



class XmippProtMultipleFSCs(ProtAnalysis3D):
    """
    Compute the FSCs between a reference volume and a set of input volumes.
    A mask can be provided and the volumes are aligned by default.
    """
    _label = 'multiple fscs'
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addParam('referenceVolume', PointerParam, label="Reference volume",
                      pointerClass='Volume',
                      help='The rest of volumes will be compared to this one')
        form.addSection(label='Input')
        form.addParam('inputVolumes', MultiPointerParam, pointerClass='Volume',
                      label="Volumes to compare",
                      help='Set of volumes to compare to the reference volume')
        form.addParam('mask', PointerParam, pointerClass='VolumeMask', allowsNull=True,
                      label="Mask",
                      help='A mask may be provided and it is applied before '
                           'comparing the different volumes')
        form.addParam('doAlign', BooleanParam, default=True,
                      label="Align volumes?",
                      help="Align volumes to reference before comparing. A local "
                           "alignment is performed so the initial orientation "
                           "of the volumes should be relatively similar")
        form.addParallelSection(threads=8, mpi=1)

#--------------------------- INSERT steps functions --------------------------------------------  
                                
    def _insertAllSteps(self):
        stepId = self._insertFunctionStep('prepareReferenceStep',
                                          self.referenceVolume.get().getObjId())
        i = 0
        for vol in self.inputVolumes:
            fnVol = getImageLocation(vol.get())
            self._insertFunctionStep('compareVolumeStep',fnVol, i, prerequisites=[stepId])
            i+=1

    def prepareReferenceStep(self,volId):
        if self.mask.hasValue():
            img=ImageHandler()
            fnMask = self._getExtraPath("mask.vol")
            img.convert(self.mask.get(), fnMask)

            referenceXdim = self.referenceVolume.get().getDim()[0]
            maskXdim= self.mask.get().getDim()[0]
            if referenceXdim!=maskXdim:
                self.runJob('xmipp_image_resize',"-i %s --dim %d"%(fnMask,referenceXdim))
            fnRef = self._getExtraPath("reference.vol")
            img.convert(self.referenceVolume.get(), fnRef)
            self.runJob("xmipp_image_operate","-i %s --mult %s"%(fnRef, fnMask))
    
    def compareVolumeStep(self,fnVol,i):
        img=ImageHandler()
        referenceXdim = self.referenceVolume.get().getDim()[0]
        fnRef = self._getExtraPath("reference.vol")
        if not os.path.exists(fnRef):
            fnRef = getImageLocation(self.referenceVolume.get())
        
        volumeXdim, _, _, _ =img.getDimensions((1,fnVol))
        fnRoot = self._getExtraPath("volume_%02d"%i)
        if referenceXdim!=volumeXdim or self.doAlign.get() or self.mask.hasValue():
            fnExtraVol = fnRoot + ".vol"
            self.runJob("xmipp_image_convert","-i %s -o %s"%(fnVol,fnExtraVol))
            if referenceXdim!=volumeXdim:
                self.runJob('xmipp_image_resize',"-i %s --dim %d"%(fnExtraVol,referenceXdim))
            if self.doAlign:
                self.runJob('xmipp_volume_align',"--i1 %s --i2 %s --apply --local"%(fnRef,fnExtraVol))
            if self.mask.hasValue():
                self.runJob("xmipp_image_operate","-i %s --mult %s"%(fnExtraVol, self._getExtraPath("mask.vol")))
            fnVol = fnExtraVol
        
        self.runJob("xmipp_resolution_fsc","--ref %s -i %s -o %s --sampling_rate %f"%\
                    (fnRef, fnVol, fnRoot+".fsc", self.referenceVolume.get().getSamplingRate()))
