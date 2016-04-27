# **************************************************************************
# *
# * Authors:         Josue Gomez Blanco (jgomez@cnb.csic.es)
# *
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

from pyworkflow.protocol.params import PointerParam
import pyworkflow.em.metadata as md

from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.em.packages.xmipp3.convert import (writeSetOfParticles,
                                                   getImageLocation)


class XmippProtSubtractProjection(ProtAnalysis3D):
    """    
    Extract the information contained in a volume to the experimental
    particles. The particles must have projection alignment in order to
    properly generate volume projections to extract the information.
    A typical case of use, is in the deletion of the capsid to the
    experimental image to only refine the genetic material.
    """
    _label = 'subtract projection'
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles", important=True,
                      help='Select the experimental particles.')
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input volume",
                      help='Select the input volume. Is desirable that the '
                           'volume was generated with the input particles.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
              label='Reference mask (optional)', allowsNull=True,
              help="The volume will be masked once the volume has been "
                   "applied the CTF of the particles.")
        
        form.addParallelSection(threads=0, mpi=0)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        imgSet = self.inputParticles.get()
        partSetId = imgSet.getObjId()
        volName = getImageLocation(self.inputVolume.get())
        self._insertFunctionStep('convertInputStep', partSetId)
        if self.refMask.get() is not None:
            self._insertFunctionStep('applyMaskStep', volName)
        self._insertFunctionStep('volumeProjectStep', volName)
        self._insertFunctionStep('removeStep')
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions -------------------------------
    def convertInputStep(self, particlesId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        imgSet = self.inputParticles.get()
        writeSetOfParticles(imgSet, self._getInputParticlesFn())
    
    def applyMaskStep(self, volName):
        params = ' -i %s --mult %s -o %s' % (volName,
                                             self.refMask.get().getFileName(),
                                             self._getOutputMap())
        self.runJob('xmipp_image_operate', params)
    
    def volumeProjectStep(self, volName):
        from convert import createParamPhantomFile
        
        imgSet = self.inputParticles.get()
        phantomFn = self._getExtraPath('params')
        pathParticles = self._getInputParticlesFn()
        dimX, _, _ = imgSet.getDim()
        
        createParamPhantomFile(phantomFn, dimX, pathParticles,
                               imgSet.isPhaseFlipped(), False)
        
        if self.refMask.get() is not None:
            volumeFn = self._getOutputMap()
        else:
            volumeFn = volName
        
        params =  ' -i %s' % volumeFn
        params += ' --params %s' % phantomFn
        params += ' -o %s' % self._getOutputRefsFn()
        params += ' --sampling_rate % 0.3f' % imgSet.getSamplingRate()
        params += ' --method fourier 2 0.5'
        self.runJob('xmipp_phantom_project', params)
    
    def removeStep(self):
        self._removeAlignLabels(self._getInputParticlesFn())
        self._removeAlignLabels(self._getOutputRefsFn())
        
        params = ' -i %s --minusAdjusted %s -o %s' % (self._getInputParticlesFn(),
                                              self._getOutputRefsFn(),
                                              self._getFinalParts())
        self.runJob('xmipp_image_operate', params)
    
    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        outImgSet = self._createSetOfParticles()
        outImgsFn = self._getFinalMDParts()
         
        outImgSet.copyInfo(imgSet)
        outImgSet.setAlignmentProj()
        outImgSet.copyItems(imgSet,
                            updateItemCallback=self._updateItem,
                            itemDataIterator=md.iterRows(outImgsFn))
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles, outImgSet)
    
    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        imgSet = self.inputParticles.get()
        vol = self.inputVolume.get()
        xImg = imgSet.getXDim()
        xVol = vol.getXDim()
        if not imgSet.getSamplingRate() == vol.getSamplingRate():
            validateMsgs.append("The sampling rate of your SetOfParticles"
                                " and your volume *MUST be equal* ")
        if not xImg == xVol:
            validateMsgs.append(" The dimensions of your particles and "
                                " your volume *MUST BE EQUAL*")
        return validateMsgs
    
    def _summary(self):
        summary = ["Input particles:  %s" % self.inputParticles.get().getNameId()]
        summary.append("-----------------")
        return summary
    
    def _methods(self):
        messages = []
        return messages
    
    #--------------------------- UTILS functions -------------------------------
    def _getInputParticlesFn(self):
        return self._getExtraPath('input_particles.xmd')
    
    def _getOutputRefsFn(self):
        return self._getExtraPath('reference_particles.xmd')
    
    def _getFinalParts(self):
        return self._getExtraPath('particles.stk')
    
    def _getFinalMDParts(self):
        return self._getExtraPath('particles.xmd')
    
    def _getOutputMap(self):
        return self._getExtraPath('masked_map.vol')
    
    def _removeAlignLabels(self, filename):
        mData = md.MetaData(filename)
        mData.removeLabel(md.MDL_ANGLE_ROT)
        mData.removeLabel(md.MDL_ANGLE_TILT)
        mData.removeLabel(md.MDL_ANGLE_PSI)
        mData.removeLabel(md.MDL_SHIFT_X)
        mData.removeLabel(md.MDL_SHIFT_Y)
        mData.removeLabel(md.MDL_SHIFT_Z)
        mData.write(filename)
    
    def _updateItem(self, item, row):
        from convert import xmippToLocation
        newFn = row.getValue(md.MDL_IMAGE)
        newLoc = xmippToLocation(newFn)
        item.setLocation(newLoc)
