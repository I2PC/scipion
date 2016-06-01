# **************************************************************************
# *
# * Authors:         Josue Gomez Blanco (jgomez@cnb.csic.es)
# *                  Roberto Marabini   (roberto@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam

from pyworkflow.em.protocol import ProtOperateParticles
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.em.packages.xmipp3.convert import (XmippMdRow, particleToRow,
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

    CORRECT_NONE = 0
    CORRECT_FULL_CTF = 1
    CORRECT_PHASE_FLIP = 2

    def __init__(self, *args, **kwargs):
        ProtOperateParticles.__init__(self, *args, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

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
        form.addParam('projType', EnumParam, default=self.CORRECT_NONE,
                      choices=['NO','full CTF','abs(CTF) (phase flip)'],
                      label='apply CTF?',
                      help='apply CTF to the reference gallery of projections')
        form.addParam('normalize', BooleanParam, default=False,
                   label='normalize', help = "adjust grey scale range so experimental data and synthetic data match")

        form.addParallelSection(threads=3, mpi=0)

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        #divide work in several threads
        numberOfThreads = self.numberOfThreads.get()-1
        if numberOfThreads == 0:
            numberOfThreads = 1
        partSet = self.inputParticles.get()
        nPart = len(partSet)
        numberOfThreads = min(numberOfThreads, nPart)
        samplingRate = partSet.getSamplingRate()
        inputVolume = getImageLocation(self.inputVolume.get())
        mdFn = self._getInputParticlesFn()
        convertId = self._insertFunctionStep('convertInputStep',
                                            inputVolume)

        deps = [convertId]

        # Create xmipp metadata
        # partSet
        #writeSetOfParticles(partSet, mdFn, blockName="images")
        groupSize = nPart/numberOfThreads
        groupRemainder = nPart%numberOfThreads

        depsOutPut=[]
        for thread in range(0, numberOfThreads):
            start = long(thread * groupSize+1)
            end = long(thread * groupSize+groupSize)
            if thread == (numberOfThreads-1):
                end += groupRemainder
            idStep = self._insertFunctionStep('projectStep', start, end,
                                              samplingRate, thread,
                                              prerequisites=deps)
            depsOutPut.append(idStep)
        self._insertFunctionStep('createOutputStep', prerequisites=depsOutPut)

    #--------------------------- STEPS functions -------------------------------
    def convertInputStep(self, volName):
        # Read volume and convert to DOUBLE
        self.vol = xmipp.Image(volName)
        self.vol.convert2DataType(xmipp.DT_DOUBLE)
        # Mask volume if needed
        if self.refMask.get() is not None:
            maskName = getImageLocation(self.refMask.get())
            self.mask = xmipp.Image(maskName)
            self.mask.convert2DataType(xmipp.DT_DOUBLE)
            self.vol.inplaceMultiply(self.mask)
        padding = 2
        maxFreq = 0.5
        splineDegree = 2
        ###
        self.fourierProjectVol = xmipp.FourierProjector(self.vol, padding, maxFreq, splineDegree)
        ###
        partSet = self.inputParticles.get()
        nPart = len(partSet)
        numberOfTasks = min(nPart, max(self.numberOfThreads.get()-1, 1))
        taskSize = nPart / numberOfTasks
        md = xmipp.MetaData()

        # Convert angles and shifts from volume system of coordinates to
        # projection system of coordinates
        mdCount = 0

        for index, part in enumerate(partSet):
            objId = md.addObject()
            imgRow = XmippMdRow()
            particleToRow(part, imgRow)
            shifts, angles = geometryFromMatrix(part.getTransform().getMatrix(), True)

            imgRow.setValue(xmipp.MDL_SHIFT_X,-shifts[0])
            imgRow.setValue(xmipp.MDL_SHIFT_Y,-shifts[1])
            imgRow.setValue(xmipp.MDL_SHIFT_Z,0.)
            imgRow.setValue(xmipp.MDL_ANGLE_ROT,angles[0])
            imgRow.setValue(xmipp.MDL_ANGLE_TILT,angles[1])
            imgRow.setValue(xmipp.MDL_ANGLE_PSI,angles[2])

            imgRow.writeToMd(md, objId)

            # Write a new metadata every taskSize number of elements
            # except in the last chunk where we want to add also the
            # remainder and the condition is the last element
            if ((index % taskSize == taskSize-1 and
                 mdCount < numberOfTasks-1) or index == nPart-1):
                md.write(self._getInputParticlesSubsetFn(mdCount))
                md.clear()
                mdCount += 1

        x, y, _ = partSet.getDim()
        xmipp.createEmptyFile(self._getProjGalleryFn(), x, y, 1, nPart)

    def projectStep(self, start, end, samplingRate, threadNumber):
        # Project
        md = xmipp.MetaData(self._getInputParticlesSubsetFn(threadNumber))
        ##
        projection = xmipp.Image()
        projection.setDataType(xmipp.DT_DOUBLE)
        ##
        for id in md:
            rot  = md.getValue(xmipp.MDL_ANGLE_ROT,  id)
            tilt = md.getValue(xmipp.MDL_ANGLE_TILT, id)
            psi  = md.getValue(xmipp.MDL_ANGLE_PSI,  id)

            ##projection =self.vol.projectVolumeDouble(rot, tilt, psi)
            self.fourierProjectVol.projectVolume(projection, rot, tilt, psi)
            ##
            # Apply CTF
            if self.projType == self.CORRECT_NONE:
                pass
            elif self.projType == self.CORRECT_FULL_CTF:
                projection.applyCTF(md, samplingRate, id, False)
            elif self.projType == self.CORRECT_PHASE_FLIP:
                projection.applyCTF(md, samplingRate, id, True)
            else:
                raise Exception("ERROR: Unknown projection mode: %d" % self.projType)

            # Shift image
            projection.applyGeo(md,id,True,False)#onlyapplyshist, wrap
            ih = ImageHandler()
            expProj = ih.read(md.getValue(xmipp.MDL_IMAGE, id))
            expProj.convert2DataType(xmipp.DT_DOUBLE)
            # Subtract from experimental and write result
            projection.resetOrigin()
            if self.normalize:
                expProj = expProj.adjustAndSubtract(projection)
            else:
                expProj.inplaceSubtract(projection)

            expProj.write( self._getProjGalleryIndexFn(id+start-1))

    def createOutputStep(self):
        partSet = self.inputParticles.get()
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(partSet)
        outImgSet.setAlignmentProj()
        for index, part in enumerate(partSet):
            part2 = part.clone()
            part2.setLocation(index+1, self._getProjGalleryFn())
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
        #if imgSet.isPhaseFlipped() and self.projType != self.CORRECT_NONE:
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

    def _getInputParticlesSubsetFn(self, setNumber):
        return self._getExtraPath('input_particles_thr_%05d.xmd'%setNumber)

    def _getProjGalleryFn(self):
        return self._getExtraPath('projGallery.stk')

    def _getProjGalleryIndexFn(self,index):
        return "%d@%s"%(index, self._getProjGalleryFn())

