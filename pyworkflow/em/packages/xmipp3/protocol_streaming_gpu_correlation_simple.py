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

from pyworkflow.em import SetOfParticles, SetOfClasses2D, ALIGN_2D, ALIGN_NONE
from pyworkflow.em.protocol import ProtAlign2D
import pyworkflow.em.metadata as md
import pyworkflow.protocol.params as params
from pyworkflow.em.metadata.utils import iterRows, getSize
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles, \
    xmippToLocation, rowToAlignment, writeSetOfClasses2D, particleToRow, \
    rowToParticle
from shutil import copy
from os.path import join, exists, getmtime
from os import mkdir, remove
from datetime import datetime
from pyworkflow.utils import prettyTime
from pyworkflow.object import Set
from pyworkflow.protocol.constants import STATUS_NEW
from xmipp3 import XmippMdRow
import time
from pyworkflow.em.data import Particle, Class2D



REF_CLASSES = 0
REF_AVERAGES = 1

class HashTableDict:
    def __init__(self, Ndict=100):
        self.Ndict = Ndict
        self.dict = [{}]*Ndict

    def isItemPresent(self, idx):
        return idx in self.dict[idx % self.Ndict]

    def pushItem(self, idx):
        idxDict = idx % self.Ndict
        if not idx in self.dict[idxDict]:
            self.dict[idxDict][idx]=1



class XmippProtStrGpuCrrSimple(ProtAlign2D):
    """ Aligns a set of particles in streaming using the GPU Correlation algorithm. """
    _label = 'align with GPU Correlation in streaming'


    # --------------------------- DEFINE param functions -----------------------
    def _defineAlignParams(self, form):
        form.addParam('useAsRef', params.EnumParam,
                      choices=['Classes', 'Averages'],
                      default=0, important=True,
                      label='Set to use as reference images',
                      display=params.EnumParam.DISPLAY_LIST,
                      help="Select which kind of set to use as reference images "
                           "for the classification of the new particles.")
        form.addParam('inputClasses', params.PointerParam,
                      pointerClass='SetOfClasses2D',
                      condition='useAsRef==%d' % REF_CLASSES,
                      label="Set of reference classes",
                      help='Set of classes that will serve as reference for '
                           'the classification')
        form.addParam('inputAverages', params.PointerParam,
                      pointerClass='SetOfAverages',
                      condition='useAsRef==%d' % REF_AVERAGES,
                      label="Set of reference averages",
                      help='Set of averages that will serve as reference for '
                           'the classification')
        form.addParam('maximumShift', params.IntParam, default=10,
                      label='Maximum shift (px):')
        form.addParam('keepBest', params.IntParam, default=2,
                      label='Number of best images:',
                      help='Number of the best images to keep for every class')
        form.addParallelSection(threads=0, mpi=0)


    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        """" Insert the steps to call cuda correlation program"""

        self.listInFn = []
        self.listOutFn = []
        self.doneListFn = []
        self.imgsRef = self._getExtraPath('imagesRef.xmd')
        self.htAlreadyProcessed = HashTableDict()

        self.last_time = time.time()

        self._loadInputList()
        if self.useAsRef.get() == 0:
            classId = self.inputClasses.get()
        else:
            classId = self.inputAverages.get()


        deps = []
        self._insertFunctionStep('convertAveragesStep', classId)
        deps = self._insertStepsForParticles(deps)

        self._insertFunctionStep('createOutputStep',
                                 prerequisites=deps, wait=True)

    def _insertStepsForParticles(self, deps):
        stepIdClassify = self._insertFunctionStep('classifyStep',
                                                  prerequisites= deps)
        deps.append(stepIdClassify)
        return deps

    # --------------------------- STEPS functions ------------------------------
    def convertAveragesStep(self, classId):

        if self.useAsRef==REF_CLASSES:
            setOfClasses = self.inputClasses
        else:
            setOfClasses = self.inputAverages

        if self.useAsRef == REF_CLASSES:
            writeSetOfClasses2D(setOfClasses.get(), self.imgsRef,
                                writeParticles=True)
        else:
            writeSetOfParticles(setOfClasses.get(), self.imgsRef)

    def classifyStep(self):
        inputImgs = self._getInputFn()
        writeSetOfParticles(self.listOfParticles, inputImgs,
                            alignType=ALIGN_NONE)

        for p in self.listOfParticles:
            partId = p.getObjId()
            self.htAlreadyProcessed.pushItem(partId)

        # Calling program xmipp_cuda_correlation
        outImgs, clasesOut = self._getOutputsFn()
        self._params = {'imgsRef': self.imgsRef,
                        'imgsExp': inputImgs,
                        'outputFile': outImgs,
                        'keepBest': self.keepBest.get(),
                        'maxshift': self.maximumShift.get(),
                        'outputClassesFile': clasesOut,
                        }

        args = ('-i_ref %(imgsRef)s -i_exp %(imgsExp)s -o %(outputFile)s '
                '--keep_best %(keepBest)d --maxShift %(maxshift)d '
                '--simplifiedMd --classify %(outputClassesFile)s')
        self.runJob("xmipp_cuda_correlation", args % self._params)

    # ------ Methods for Streaming 2D Classification --------------
    def _stepsCheck(self):
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        """ Check if there are new particles to be processed and add
        the necessary steps."""
        particlesFile = self.inputParticles.get().getFileName()

        now = datetime.now()
        self.lastCheck = getattr(self, 'lastCheck', now)
        mTime = datetime.fromtimestamp(getmtime(particlesFile))
        self.debug('Last check: %s, modification: %s'
                   % (prettyTime(self.lastCheck),
                      prettyTime(mTime)))

        # If the input have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime and hasattr(self, 'listOfParticles'):
            return None

        self.lastCheck = now
        outputStep = self._getFirstJoinStep()

        # Open input and close it as soon as possible
        self._loadInputList()

        fDeps=[]
        fDeps = self._insertStepsForParticles(fDeps)
        if outputStep is not None:
            outputStep.addPrerequisites(*fDeps)
        self.updateSteps()

    def _checkNewOutput(self):
        """ Check for already done files and update the output set. """

        # Check for newly done items
        newDone = self._readDoneList()

        # We have finished when there is not more inputs (stream closed)
        # and the number of processed particles is equal to the number of inputs
        self.finished = (self.isStreamClosed == Set.STREAM_CLOSED
                         and len(newDone)==0)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if newDone:
            #self._updateOutputCoordSet(newDone, streamMode)
            self._updateOutputSetOfClasses(newDone, streamMode)
            self.doneListFn += newDone



        elif not self.finished:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        #self.debug('   finished: %s ' % self.finished)
        #self.debug('        self.streamClosed (%s) AND' % streamMode)
        #self.debug('        allDone (%s) == len(self.listOfMics (%s)'
        #           % (allDone, nMics))

        if self.finished:  # Unlock createOutputStep if finished all jobs
            #self._updateStreamState(streamMode)
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)






        #outSet = self._loadOutputSet(SetOfClasses2D, 'classes2D.sqlite')
        #self._updateOutputSetOfClasses('outputClasses', outSet, streamMode)
        #if (exists(self._getPath('classes2D.sqlite'))):
        #    outSet.close()





    def createOutputStep(self):
        pass

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        if self.useAsRef==REF_CLASSES:
            refImage = self.inputClasses.get()
        else:
            refImage = self.inputAverages.get()
        [x1, y1, z1] = refImage.getDimensions()
        [x2, y2, z2] = self.inputParticles.get().getDim()
        if x1 != x2 or y1 != y2 or z1 != z2:
            errors.append('The input images and the reference images '
                          'have different sizes')
        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputClasses'):
            summary.append("Output alignment not ready yet.")
        else:
            summary.append("Input Particles: %s"
                           % self.inputParticles.get().getSize())
            summary.append("Aligned with reference images: %s"
                           % self.inputClasses.get().getSize())
        return summary

    def _citations(self):
        return ['Sorzano2010a']

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputClasses'):
            methods.append("Output alignment not ready yet.")
        else:
            methods.append(
                "We aligned images %s with respect to the reference image set "
                "%s using Xmipp CUDA correlation"
                % (self.getObjectTag('inputParticles'),
                   self.getObjectTag('inputClasses')))

        return methods

    # --------------------------- UTILS functions ------------------------------
    def _loadInputList(self):
        """ Load the input set of ctfs and create a list. """
        initial_time = time.time()
        particlesSet = self._loadInputParticleSet()
        self.isStreamClosed = particlesSet.getStreamState()
        self.listOfParticles = []
        for p in particlesSet:
            idx = p.getObjId()
            if not self.htAlreadyProcessed.isItemPresent(idx):
                #self.listOfParticles.append(idx)
                # AJ new to make the list of particles instead of list of ids
                newPart = p.clone()
                self.listOfParticles.append(newPart)
        particlesSet.close()
        self.debug("Closed db.")

    def _loadInputParticleSet(self):
        partSetFn = self.inputParticles.get().getFileName()
        updatedSet = SetOfParticles(filename=partSetFn)
        copyPartSet = SetOfParticles()
        updatedSet.loadAllProperties()
        copyPartSet.copy(updatedSet)
        updatedSet.close()
        return copyPartSet

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == 'createOutputStep':
                return s
        return None


    def _readDoneList(self):
        return [fn for fn in self.listOutFn if fn not in self.doneListFn]

    def _loadOutputSet(self, SetClass, baseName):

        setFile = self._getPath(baseName)
        if exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
            outputSet.setStreamState(outputSet.STREAM_OPEN)
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        inputs = self._loadInputParticleSet()
        outputSet.copyInfo(inputs)
        if SetClass is SetOfClasses2D:
            outputSet.setImages(inputs)
        return outputSet


    def _updateOutputSetOfClasses(self, outFnDone, streamMode):
        outputName = 'outputClasses'
        outputClasses = getattr(self, outputName, None)
        firstTime = True

        if outputClasses is None:
            outputClasses = self._createSetOfClasses2D(self.inputParticles.get())
        else:
            firstTime = False
            outputClasses.enableAppend()
            # for class2d in outputClasses:
            #     class2d.enableAppend()

        self._fillClassesFromMd(outFnDone, outputClasses, firstTime)
        self._updateOutputSet(outputName, outputClasses, streamMode)

        if firstTime:
            self._defineSourceRelation(self.inputParticles, outputClasses)




        # outputSet.setStreamState(state)
        # if self.hasAttribute(outputName):
        #     print("En el if")
        #     outputSet.enableAppend()
        #     self._fillClassesFromLevel(outputSet, False)
        #     outputSet.write()  # Write to commit changes
        #     outputAttr = getattr(self, outputName)
        #     # Copy the properties to the object contained in the protocol
        #     outputAttr.copy(outputSet, copyId=False)
        #     # Persist changes
        #     self._store(outputAttr)
        # else:
        #     print("En el else")
        #     self._fillClassesFromLevel(outputSet, True)
        #     self._defineOutputs(**{outputName: outputSet})
        #     self._defineSourceRelation(self.inputParticles, outputSet)
        #     self._store(outputSet)
        #
        # # Close set database to avoid locking it
        # outputSet.close()

    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(md.MDL_REF))
        item.setTransform(rowToAlignment(row, ALIGN_2D))

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn = self._classesInfo[classId]
            item.setAlignment2D()
            rep = item.getRepresentative()
            rep.setLocation(index, fn)
            rep.setSamplingRate(self.inputParticles.get().getSamplingRate())


    def _loadClassesInfo(self):
        """ Read some information about the produced 2D classes
        from the metadata file.
        """
        self._classesInfo={}
        for row in md.iterRows(self.imgsRef):
            index, fn = xmippToLocation(row.getValue(md.MDL_IMAGE))
            self._classesInfo[index] = (index, fn)

    def _fillClassesFromMd(self, outFnDone, outputClasses, firstTime):

        for outFn in outFnDone:
            myFileClasses = outFn.replace('images', 'classes')
            self._loadClassesInfo()
            if firstTime:
                inputSet = self._loadInputParticleSet()
                mdClass = md.MetaData("classes@" + myFileClasses)
                rows = iterRows(mdClass)
                for myClass in rows:
                    ref = myClass.getValue(md.MDL_REF)
                    classItem = Class2D(objId=ref)
                    classItem.setStreamState(Set.STREAM_OPEN)
                    rep = Particle()
                    classItem.setRepresentative(rep)
                    classItem.copyInfo(inputSet)
                    classItem.setAcquisition(inputSet.getAcquisition())
                    self._updateClass(classItem)
                    outputClasses.append(classItem)
                    mdImages = md.MetaData("class%06d_images@" % ref + myFileClasses)
                    images = iterRows(mdImages)
                    for image in images:
                        if image.getValue(md.MDL_ENABLED) == 1:
                            rowImage = rowToParticle(image)
                            classItem.append(rowImage)
                    outputClasses.update(classItem)
                    firstTime = False




    def _fillClassesFromLevel(self, outputSet, firstTime):
        """ Create the SetOfClasses2D from a given iteration. """

        myFileParticles = self._getExtraPath('last_images.xmd')
        myFileClasses = self._getExtraPath('last_classes.xmd')

        self._loadClassesInfo(myFileClasses)

        if firstTime:
            #iterator = md.SetMdIterator(myFileParticles, sortByLabel=md.MDL_ITEM_ID,
            #                            updateItemCallback=self._updateParticle,
            #                            skipDisabled=True)
            ## itemDataIterator is not necessary because, the class SetMdIterator
            ## contain all the information about the metadata
            #clsSet.classifyItems(updateItemCallback=iterator.updateItem,
            #                     updateClassCallback=self._updateClass)
            inputSet = self._loadInputParticleSet()
            mdClass = md.MetaData("classes@" + myFileClasses)
            rows = iterRows(mdClass)
            for myClass in rows:
                ref = myClass.getValue(md.MDL_REF)
                classItem = Class2D(objId=ref)
                classItem.setStreamState(Set.STREAM_OPEN)
                rep = Particle()
                classItem.setRepresentative(rep)
                classItem.copyInfo(inputSet)
                classItem.setAcquisition(inputSet.getAcquisition())
                self._updateClass(classItem)
                outputSet.append(classItem)
                mdImages = md.MetaData("class%06d_images@"%ref + myFileClasses)
                images = iterRows(mdImages)
                for image in images:
                    if image.getValue(md.MDL_ENABLED)==1:
                        rowImage = rowToParticle(image)
                        classItem.append(rowImage)
                outputSet.update(classItem)

        else:
            for myClass in outputSet.iterItems():
                myClass.enableAppend()
                ref = myClass.getObjId()
                print("ref",ref)
                mdImages = md.MetaData("class%06d_images@" % ref + myFileClasses)
                images = iterRows(mdImages)
                for image in images:
                    if image.getValue(md.MDL_ENABLED) == 1:
                        rowImage = rowToParticle(image)
                        myClass.append(rowImage)
                outputSet.update(myClass)
                #outputSet._insertItem(myClass)

    def _getUniqueFn(self, basename, list):
        if list == []:
            fn = basename + "_1.xmd"
        else:
            number = int(list[-1].split("_")[-1].split(".")[0]) + 1
            fn = basename + "_%s.xmd" % number
        list.append(fn)
        return fn

    def _getInputFn(self):
        basename = self._getExtraPath('imagesExp')
        return self._getUniqueFn(basename, self.listInFn)

    def _getOutputsFn(self):
        nameImages = self._getExtraPath('general_images')
        imagesFn = self._getUniqueFn(nameImages, self.listOutFn)
        classesFn = imagesFn.replace('images', 'classes')
        return imagesFn, classesFn
