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

from pyworkflow.protocol.constants import STEPS_PARALLEL
import pyworkflow.protocol.params as params

import pyworkflow.em as em
import pyworkflow.em.metadata as md



class XmippProtMultipleFSCs(em.ProtAnalysis3D):
    """
    Compute the FSCs between a reference volume and a set of input volumes.
    A mask can be provided and the volumes are aligned by default.
    """
    _label = 'multiple fscs'

    def __init__(self, **args):
        em.ProtAnalysis3D.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('referenceVolume', params.PointerParam,
                      pointerClass='Volume',
                      label="Reference volume",
                      help='The rest of volumes will be compared to this one')

        form.addParam('inputVolumes', params.MultiPointerParam,
                      pointerClass='Volume',
                      label="Volumes to compare",
                      help='Set of volumes to compare to the reference volume')
        form.addParam('mask', params.PointerParam,
                      pointerClass='VolumeMask', allowsNull=True,
                      label="Mask",
                      help='A mask may be provided and it is applied before '
                           'comparing the different volumes')
        form.addParam('doAlign', params.BooleanParam, default=True,
                      label="Align volumes?",
                      help="Align volumes to reference before comparing. A local "
                           "alignment is performed so the initial orientation "
                           "of the volumes should be relatively similar")
        form.addParallelSection(threads=8, mpi=1)

#--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        stepId = self._insertFunctionStep('prepareReferenceStep',
                                          self.referenceVolume.get().getObjId())
        allVols = []
        for i, vol in enumerate(self.inputVolumes):
            volId = self._insertFunctionStep('compareVolumeStep',
                                             vol.get().getLocation(), i+1,
                                             prerequisites=[stepId])
            allVols.append(volId)

        self._insertFunctionStep('createOutputStep', prerequisites=allVols)

    def _resizeVolume(self, volFn):
        """ Resize input volume if not of the same size of referenceVol """
        refDim = self.referenceVolume.get().getXDim()
        volDim = em.ImageHandler().getDimensions(volFn)
        if refDim != volDim:
            self.runJob('xmipp_image_resize',
                        "-i %s --dim %d" % (volFn, refDim))

    def _maskVolume(self, volFn):
        """ Mask input volume multiplying by mask.vol. """
        self.runJob("xmipp_image_operate",
                    "-i %s --mult %s" % (volFn, self._getExtraPath("mask.vol")))

    def prepareReferenceStep(self,volId):
        inputMask = self.mask.get()

        ih = em.ImageHandler()
        fnRef = self._getExtraPath("reference.vol")
        ih.convert(self.referenceVolume.get(), fnRef)

        if inputMask is not None:
            fnMask = self._getExtraPath("mask.vol")
            ih.convert(self.mask.get(), fnMask)
            self._resizeVolume(fnMask)
            self._maskVolume(fnRef)

    def compareVolumeStep(self, volLoc, i):
        ih = em.ImageHandler()
        fnRef = self._getExtraPath("reference.vol")
        sampling = self.referenceVolume.get().getSamplingRate()
        fnRoot = self._getExtraPath("volume_%02d" % i)
        fnVol = fnRoot + ".vol"
        ih.convert(volLoc, fnVol)

        # Resize if the volume has different size than the reference
        self._resizeVolume(fnVol)

        if self.doAlign: # Align against the reference if selected
            self.runJob('xmipp_volume_align',
                        "--i1 %s --i2 %s --apply --local" % (fnRef, fnVol))

        if self.mask.hasValue(): # Mask volume if input mask
            self._maskVolume(fnVol)

        # Finally compute the FSC
        args = "--ref %s -i %s -o %s_fsc.xmd --sampling_rate %f" % (fnRef, fnVol,
                                                                fnRoot, sampling)
        self.runJob("xmipp_resolution_fsc", args)

    def createOutputStep(self):
        fscSet = self._createSetOfFSCs()

        for i, vol in enumerate(self.inputVolumes):
            index = i + 1
            fnFsc = self._getExtraPath("volume_%02d_fsc.xmd" % index)
            mdFsc = md.MetaData(fnFsc)
            fscLabel = vol.get().getObjLabel() or 'FSC %d' % index
            fsc = em.FSC(objLabel=fscLabel)
            fsc.loadFromMd(mdFsc, md.MDL_RESOLUTION_FREQ, md.MDL_RESOLUTION_FRC)
            fscSet.append(fsc)

        self._defineOutputs(outputFSCs=fscSet)

        self._defineSourceRelation(self.referenceVolume, fscSet)
        for i, vol in enumerate(self.inputVolumes):
            self._defineSourceRelation(vol, fscSet)

