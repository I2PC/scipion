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
from pyworkflow.em.packages.xmipp3.xmipp3 import xmippExists
from pyworkflow.protocol.params import PointerParam, EnumParam
import pyworkflow.em.metadata as md

from pyworkflow.em.protocol import ProtOperateParticles
from pyworkflow.em.packages.xmipp3.convert import (writeSetOfParticles,
                                                   getImageLocation,
                                                   geometryFromMatrix)
import xmipp
from pyworkflow.em import ImageHandler

class XmippProtSubtractProjection(ProtOperateParticles):
    """    
    Subtract volume projections from the experimental particles.
    The particles must have projection alignment in order to
    properly generate volume projections.

    An example of usage is to delete the virus capsid to
    refine only the genetic material.
    """
    _label = 'subtract projection'
    _CorrectNone = 0
    _correctFullCTF = 1
    _correctPhaseFlip =2
    _projGalleryFn = "projGallery.stk"
    _metadata      = "inputParticles.xmd"

    def __init__(self, *args, **kwargs):
        ProtOperateParticles.__init__(self, *args, **kwargs)

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
        form.addParam('projType', EnumParam, default=self._CorrectNone,
                      choices=['NO','full CTF','abs(CTF) (phase flip)'],
                      label='apply CTF?',
                      help='apply CTF to the reference gallery of projections')

        form.addParallelSection(threads=3, mpi=0)

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        deps = []
        partSet = self.inputParticles.get()
        samplingRate = partSet.getSamplingRate()
        inputVolume = getImageLocation(self.inputVolume.get())
        mdFn = self._getExtraPath(self._metadata)
        volumeId = self._insertFunctionStep('initVolumeStep', inputVolume, partSet, mdFn)
        deps.append(volumeId)

        projGalleryFn = self._getExtraPath(self._projGalleryFn)
        galleryId = self._insertFunctionStep('createEmptyFileStep', projGalleryFn)
        deps.append(galleryId)

        #create xmipp metadata
        writeSetOfParticles(partSet, mdFn, blockName="images")
        self.md = xmipp.MetaData(mdFn)

        #insert one step per projection direction
        for index, part in enumerate(partSet):
            indexPart, expProjFilename = part.getLocation()
            shifts, angles = geometryFromMatrix(part.getTransform().getMatrix(), True)
            self._insertFunctionStep('projectStep', index+1, shifts, angles,
                                     samplingRate,expProjFilename, projGalleryFn,
                                     prerequisites=deps)
        self._insertFunctionStep('createOutputStep', projGalleryFn)

    #--------------------------- STEPS functions -------------------------------
    def initVolumeStep(self,volName, partSet, mdFn):
        #read volume
        self.vol = xmipp.Image(volName)
        self.vol.convert2DataType(xmipp.DT_DOUBLE)
        #self.xdim, self.ydim, self.zdim, self.n = self.vol.getDimensions()

        #mask volume if needed
        if self.refMask.get() is not None:
            maskName = getImageLocation(self.refMask.get())
            self.mask = xmipp.Image(maskName)
            self.mask.convert2DataType(xmipp.DT_DOUBLE)
            self.vol.inplaceMultiplyImg(self.mask)

    def createEmptyFileStep(self,projGalleryFn):
        n = len(self.inputParticles.get())
        x, y, z = self.inputParticles.get().getDimensions()
        xmipp.createEmptyFile(projGalleryFn,x,y,1,n)
        #createEmptyFile(PyString_AsString(input),Xdim,Ydim,Zdim,Ndim,true,WRITE_REPLACE);

    def projectStep(self,i, shifts, angles, samplingRate, expProjFilename, projGalleryFn):
        #project
        #shifts, angles = geometryFromMatrix(part.getTransform().getMatrix(), True)
        self.projection =self.vol.projectVolumeDouble(angles[0], angles[1], angles[2])
        self.md.setValue(xmipp.MDL_SHIFT_X,-shifts[0],i)
        self.md.setValue(xmipp.MDL_SHIFT_Y,-shifts[1],i)
        self.md.setValue(xmipp.MDL_SHIFT_Z,0.,i)
        #apply CTF
        if self.projType == self._CorrectNone:
            pass
        elif self.projType == self._correctFullCTF:
            self.projection.applyCTF(self.md, samplingRate, i, False)
        elif self.projType == self._correctPhaseFlip:
            self.projection.applyCTF(self.md, samplingRate, i, True)
        else:
            print ("ERROR: wrong projtype: ", self.projType)
            exit(1)
        #shift image
        self.projection.applyGeo(self.md,i,True,False)#onlyapplyshist, wrap
        #self.projection.write("%d@%s"%(i,projGalleryFn))
        ih = ImageHandler()
        print("expProjFilename", "%d@%s"%(i,expProjFilename))
        expProj = ih.read("%d@%s"%(i,expProjFilename))
        expProj.convert2DataType(xmipp.DT_DOUBLE)
        #substract from experimental
        #write result
        print expProj.getDimensions()
        print self.projection.getDimensions()
        self.projection.resetOrigin()
        expProj.inplaceSubtract(self.projection) #0, -64
        expProj.write("%d@%s"%(i,projGalleryFn))

        #do not forget to handle shifts
        #applyTransforMatScipion

    def DELETEconvertInputStep(self, particlesId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        imgSet = self.inputParticles.get()
        writeSetOfParticles(imgSet, self._getInputParticlesFn())

    def DELETEapplyMaskStep(self, volName, maskName):
        vol = xmipp.Image()
        params = ' -i %s --mult %s -o %s' % (volName,
                                             maskName,
                                             self._getOutputMap())
        self.runJob('xmipp_image_operate', params)

    def DELETEvolumeProjectStep(self, volName):
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

    def DELETEremoveStep(self):
        self._removeAlignLabels(self._getInputParticlesFn())
        self._removeAlignLabels(self._getOutputRefsFn())

        params = ' -i %s --minusAdjusted %s -o %s' % (self._getInputParticlesFn(),
                                                      self._getOutputRefsFn(),
                                                      self._getFinalParts())
        self.runJob('xmipp_image_operate', params)

    def createOutputStep(self, projGalleryFn):
        partSet = self.inputParticles.get()
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(partSet)
        outImgSet.setAlignmentProj()
        for index, part in enumerate(partSet):
            #indexPart, expProjFilename = part.getLocation()
            part2 = part.clone()
            part2.setFileName(projGalleryFn)
            part2.setIndex(index+1)
            outImgSet.append(part2)

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
        #if imgSet.isPhaseFlipped() and self.projType != self._CorrectNone:
        #   validateMsgs.append("your images and volume are phase flipped therefore does not make sense to apply abs(CTF)")

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
