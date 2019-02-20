# ******************************************************************************
# *
# * Authors:     Amaya Jimenez Moreno (ajimenez@cnb.csic.es)
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
# ******************************************************************************

from pyworkflow.em.protocol import ProtClassify2D
import pyworkflow.em.metadata as md
import pyworkflow.protocol.params as params
from pyworkflow.em.data import Transform
from xmipp import MD_APPEND
from pyworkflow.em.packages.xmipp3.convert import rowToAlignment, \
    alignmentToRow, rowToParticle
from convert import writeSetOfClasses2D
from pyworkflow.em import ALIGN_2D
from pyworkflow.em.data import Class2D, Particle, Coordinate
import numpy as np


class XmippProtReAlignClasses(ProtClassify2D):
    """ Realignment of un-centered classes. """
    _label = 'realignment classes'

    # --------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputClasses', params.PointerParam,
                      pointerClass='SetOfClasses2D',
                      important=True,
                      label="Input Classes",
                      help='Set of classes to be realing')
        form.addParam('inputMics', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      important=True,
                      label="Set of micrographs",
                      help='Set of micrographs related to the selected input '
                           'classes')
        form.addParallelSection(threads=0, mpi=0)

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        """Mainly prepare the command line for call cuda corrrelation program"""
        self._insertFunctionStep('realignStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def realignStep(self):

        inputMdName = self._getExtraPath('inputClasses.xmd')
        writeSetOfClasses2D(self.inputClasses.get(), inputMdName,
                            writeParticles=True)

        centeredStackName = self._getExtraPath('centeredStack.stk')
        self._params = {'input': inputMdName,
                        'output': centeredStackName}
        args = ('-i %(input)s -o %(output)s --save_metadata_transform')
        self.runJob("xmipp_transform_center_image", args % self._params,
                    numberOfMpi=1)

        centeredMdName = centeredStackName.replace('stk', 'xmd')
        centeredMd = md.MetaData(centeredMdName)
        centeredStack = md.MetaData(centeredStackName)

        listName = []
        listTransform=[]
        for rowStk in md.iterRows(centeredStack):
            listName.append(rowStk.getValue(md.MDL_IMAGE))
        for rowMd in md.iterRows(centeredMd):
            listTransform.append(rowToAlignment(rowMd, ALIGN_2D))

        mdNewClasses = md.MetaData()
        for i, row in enumerate(md.iterRows(inputMdName)):
            newRow = md.Row()
            newRow.setValue(md.MDL_IMAGE, listName[i])
            refNum = row.getValue(md.MDL_REF)
            newRow.setValue(md.MDL_REF, refNum)
            classCount = row.getValue(md.MDL_CLASS_COUNT)
            newRow.setValue(md.MDL_CLASS_COUNT, classCount)
            newRow.addToMd(mdNewClasses)
        mdNewClasses.write('classes@' + self._getExtraPath('final_classes.xmd'),
                           MD_APPEND)

        mdImages = md.MetaData()
        i=0
        mdBlocks = md.getBlocksInMetaDataFile(inputMdName)
        resultMat = Transform()
        for block in mdBlocks:
            if block.startswith('class00'):
                newMat = listTransform[i]
                newMatrix = newMat.getMatrix()
                mdClass = md.MetaData(block + "@" + inputMdName)
                mdNewClass = md.MetaData()
                i+=1
                for rowIn in md.iterRows(mdClass):
                    #To create the transformation matrix (and its parameters)
                    #  for the realigned particles
                    if rowIn.getValue(md.MDL_ANGLE_PSI)!=0:
                        flag_psi=True
                    if rowIn.getValue(md.MDL_ANGLE_ROT)!=0:
                        flag_psi=False
                    inMat = rowToAlignment(rowIn, ALIGN_2D)
                    inMatrix = inMat.getMatrix()
                    resultMatrix = np.dot(newMatrix,inMatrix)
                    resultMat.setMatrix(resultMatrix)
                    rowOut=md.Row()
                    rowOut.copyFromRow(rowIn)
                    alignmentToRow(resultMat, rowOut, ALIGN_2D)
                    if flag_psi==False:
                        newAngle = rowOut.getValue(md.MDL_ANGLE_PSI)
                        rowOut.setValue(md.MDL_ANGLE_PSI, 0.)
                        rowOut.setValue(md.MDL_ANGLE_ROT, newAngle)

                    #To create the new coordinates for the realigned particles
                    inPoint = np.array([[0.],[0.],[0.],[1.]])
                    invResultMat = np.linalg.inv(resultMatrix)
                    centerPoint = np.dot(invResultMat,inPoint)
                    rowOut.setValue(md.MDL_XCOOR, rowOut.getValue(
                        md.MDL_XCOOR)+int(centerPoint[0]))
                    rowOut.setValue(md.MDL_YCOOR, rowOut.getValue(
                        md.MDL_YCOOR)+int(centerPoint[1]))

                    rowOut.addToMd(mdNewClass)
                mdNewClass.write(block + "@" + self._getExtraPath(
                    'final_classes.xmd'), MD_APPEND)
                mdImages.unionAll(mdNewClass)
        mdImages.write(self._getExtraPath('final_images.xmd'))


    def createOutputStep(self):
        inputParticles = self.inputClasses.get().getImages()
        outputClasses = self._createSetOfClasses2D(inputParticles)
        self._fillClasses(outputClasses)
        result = {'outputClasses': outputClasses}
        self._defineOutputs(**result)
        self._defineSourceRelation(self.inputClasses, outputClasses)

        outputParticles = self._createSetOfParticles()
        outputParticles.copyInfo(inputParticles)
        outputCoords = self._createSetOfCoordinates(self.inputMics.get())

        self._fillParticles(outputParticles, outputCoords, inputParticles)
        result = {'outputParticles': outputParticles}
        self._defineOutputs(**result)
        self._defineSourceRelation(self.inputClasses, outputParticles)

        result = {'outputCoordinates': outputCoords}
        self._defineOutputs(**result)
        self._defineSourceRelation(self.inputClasses, outputCoords)


    # --------------------------- UTILS functions ------------------------------

    def _fillClasses(self, outputClasses):
        """ Create the SetOfClasses2D """
        inputSet = self.inputClasses.get().getImages()
        myRep = md.MetaData('classes@' + self._getExtraPath(
            'final_classes.xmd'))

        for row in md.iterRows(myRep):
            fn = row.getValue(md.MDL_IMAGE)
            rep = Particle(fn)
            repId = row.getObjId()
            newClass = Class2D(objId=repId)
            newClass.setAlignment2D()
            newClass.copyInfo(inputSet)
            newClass.setAcquisition(inputSet.getAcquisition())
            newClass.setRepresentative(rep)
            outputClasses.append(newClass)

        i=1
        mdBlocks = md.getBlocksInMetaDataFile(self._getExtraPath(
            'final_classes.xmd'))
        for block in mdBlocks:
            if block.startswith('class00'):
                mdClass = md.MetaData(block + "@" + self._getExtraPath(
                                      'final_classes.xmd'))
                imgClassId = i
                newClass = outputClasses[imgClassId]
                newClass.enableAppend()
                for row in md.iterRows(mdClass):
                    part = rowToParticle(row)
                    newClass.append(part)
                i+=1
                newClass.setAlignment2D()
                outputClasses.update(newClass)

    def _fillParticles(self, outputParticles, outputCoords, inputParticles):
        """ Create the SetOfParticles and SetOfCoordinates"""
        myParticles = md.MetaData(self._getExtraPath('final_images.xmd'))
        outputParticles.enableAppend()
        outputCoords.enableAppend()

        #Calculating the scale that relates the coordinates with the actual
        # position in the mic
        scale = inputParticles.getSamplingRate() / \
                self.inputMics.get().getSamplingRate()
        #Dictionary with the name and id of the inpt mics
        micDictname = {}
        micDictId = {}
        for mic in self.inputMics.get():
            micKey = mic.getMicName()
            micDictname[micKey] = mic.clone()
            micKey2 = mic.getObjId()
            micDictId[micKey2] = mic.clone()

        for row in md.iterRows(myParticles):
            #To create the new particle
            p = rowToParticle(row)
            outputParticles.append(p)

            #To create the new coordinate
            newCoord = Coordinate()
            coord = p.getCoordinate()
            if coord.getMicName() is not None:
                micKey = coord.getMicName()
                micDict = micDictname
            else:
                micKey = coord.getMicId()
                micDict = micDictId
            mic = micDict.get(micKey, None)
            if mic is None:
                print "Skipping particle, key %s not found" % micKey
            else:
                newCoord.copyObjId(p)
                x, y = coord.getPosition()
                newCoord.setPosition(x * scale, y * scale)
                newCoord.setMicrograph(mic)
                outputCoords.append(newCoord)

        boxSize = inputParticles.getXDim() * scale
        outputCoords.setBoxSize(boxSize)



    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    def _summary(self):
        summary = []
        summary.append("Realignment of %s classes."
                       % self.inputParticles.get().getSize())
        return summary

