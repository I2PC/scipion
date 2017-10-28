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

from pyworkflow.em import ALIGN_NONE, SetOfParticles, SetOfClasses2D, ALIGN_2D
from pyworkflow.em.protocol import ProtAlign2D
import pyworkflow.em.metadata as md
import pyworkflow.protocol.params as params
from pyworkflow.em.metadata.utils import iterRows, getSize
import xmipp
from xmipp import MD_APPEND
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles, xmippToLocation, rowToAlignment
from shutil import copy
from os.path import join, exists, getmtime
from os import mkdir, remove, listdir
from datetime import datetime
from pyworkflow.utils import prettyTime, cleanPath
from pyworkflow.object import Set
from pyworkflow.protocol.constants import STATUS_NEW
import sys
from random import randint


class XmippProtStrGpuCrrCL2D(ProtAlign2D):
    """ Aligns a set of particles using the GPU Correlation algorithm. """
    _label = 'align with GPU Correlation in streaming'


    # --------------------------- DEFINE param functions -----------------------
    def _defineAlignParams(self, form):
        form.addParam('useReferenceImages', params.BooleanParam, default=False,
                      label='Use a Set of Reference Images ?',
                      help='If you set to *Yes*, you should provide a '
                           'set of reference images.\n'
                           'If *No*, the default generation is done by '
                           'averaging subsets of the input images.')
        form.addParam('referenceImages', params.PointerParam,
                      condition='useReferenceImages',
                      pointerClass='SetOfParticles', allowsNull=True,
                      label="Reference images",
                      help='Set of images that will serve as class reference')
        form.addParam('maximumShift', params.IntParam, default=10,
                      label='Maximum shift (px):')
        form.addParam('keepBest', params.IntParam, default=2,
                      label='Number of best images:',
                      help='Number of the best images to keep for every class')
        form.addParam('numberOfSplitIterations', params.IntParam, default=2,
                      label='Number of iterations in split stage:',
                      help='Maximum number of iterations in split stage')
        form.addParam('numberOfClassifyIterations', params.IntParam, default=4,
                      label='Number of iterations in classify stage:',
                      help='Maximum number of iterations when the classification of the whole image set is carried out')
        #form.addParam('numberOfClasses', params.IntParam, default=5,
        #              label='Number of classes:',
        #              help='Number of classes (or references) to be generated.')
        form.addParam('useAttraction', params.BooleanParam, default=True,
                      label='Allow attraction ?',
                      help='If you set to *Yes*, you allow to generate classes '
                           'with low number of images associated.\n'
                           'If *No*, all the generated classes will be '
                           'balanced.')
        form.addParallelSection(threads=0, mpi=4)


    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        """" Mainly prepare the command line for calling cuda corrrelation program"""

        self.p = 0.2
        self.percentStopClassify = 5

        self.listContinueClass=[]
        self.iterReturnSplit = 0
        self.iterReturnClass = 0
        self.numRefs = []
        self.depthSplit = 0
        self.depth = 0
        self.imgsExp = self._getExtraPath('imagesExp.xmd')
        self.exp = self.imgsExp
        self.level=0
        self.countStep = 0
        self.listNumImgs = []
        self.listNameImgs = []
        self.listRefImgs = []

        self.insertedList = []
        self.insertedListFlag = {}

        self._loadInputList()
        fDeps = self._insertProcessingStep(self.insertedList, self.insertedListFlag, self.listOfParticles, True, False)

        self._insertFunctionStep('createOutputStep', prerequisities=fDeps, wait=True)


    def _insertProcessingStep(self, insertedList, insertedListFlag, inputParticles, flagSplit, reclassification):

        print("En _insertProcessingStep",self.countStep)
        flagNewParticles = False

        self.imgsExp = self._getExtraPath('imagesExp.xmd')

        for expImgsId in inputParticles:
            if expImgsId not in insertedList:
                flagNewParticles = True
                insertedList.append(expImgsId)
                insertedListFlag[expImgsId]=self.countStep

        if flagNewParticles:
            stepId = self._insertStepsForParticles(flagSplit, reclassification)

        self.countStep+=1
        return stepId

    def levelUp(self, expImgMd):
        if getSize(expImgMd) == 0:
            print("NOOOOOOOOOOOOOOOOOOOOOOOOOOO")
            return
        self.level += 1

    def _insertStepsForParticles(self, flagSplit, reclassification):

        deps = []
        inputParticles = self.inputParticles.get()
        self.convertSet(self.imgsExp, inputParticles)
        stepIdInputmd = self._insertFunctionStep('_generateInputMd', reclassification, prerequisites=[])
        deps.append(stepIdInputmd)

        expImgMd = self._getExtraPath('inputImagesExp.xmd')
        if flagSplit:
            stepIdSplit = self._insertFunctionStep('_splitStep', expImgMd, 0, False, prerequisites=[stepIdInputmd])
            deps.append(stepIdSplit)
            stepIdLevelUp = self._insertFunctionStep('levelUp', expImgMd, prerequisites=[stepIdSplit])
            stepIdClassify = self._insertFunctionStep('_classifyStep', expImgMd, 0, prerequisites=[stepIdLevelUp])
        else:
            stepIdClassify = self._insertFunctionStep('_classifyStep', expImgMd, 0, prerequisites=[stepIdInputmd])

        stepIdCheckSplit = self._insertFunctionStep('checkSplit', expImgMd, prerequisites=[stepIdClassify])
        stepIdLevelUp = self._insertFunctionStep('levelUp', expImgMd, prerequisites=[stepIdCheckSplit])

        deps.append(stepIdLevelUp)

        return deps


    # --------------------------- STEPS functions --------------------------

    def _stepsCheck(self):
        
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        """ Check if there are new particles to be processed and add the necessary
        steps."""
        print("En _checkNewInput")
        sys.stdout.flush()

        particlesFile = self.inputParticles.get().getFileName()

        now = datetime.now()
        self.lastCheck = getattr(self, 'lastCheck', now)
        mTime = datetime.fromtimestamp(getmtime(particlesFile))
        self.debug('Last check: %s, modification: %s'
                   % (prettyTime(self.lastCheck),
                      prettyTime(mTime)))

        # Open input and close it as soon as possible
        self._loadInputList()
        #self.isStreamClosed = self.inputParticles.get().getStreamState()
        #self.listOfParticles = self.inputParticles.get()
        # If the input have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime and hasattr(self, 'listOfParticles'):
            return None

        self.lastCheck = now
        #newParticles = any(particle.getObjId() not in self.insertedList
        #             for particle in self.listOfParticles)
        newParticles = []
        for particleId in self.listOfParticles:
            if particleId not in self.insertedList:
                newParticles.append(particleId)
        outputStep = self._getFirstJoinStep()

        if len(newParticles)>0:
            print("hay new particle", newParticles)
            sys.stdout.flush()
            #print("self.insertedList",self.insertedList)
            #sys.stdout.flush()
            print("self.listOfParticles", len(self.listOfParticles))
            sys.stdout.flush()
            reclassification = True
            fDeps = self._insertProcessingStep(self.insertedList, self.insertedListFlag, newParticles, False, reclassification)
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()


    def _checkNewOutput(self):
        """ Check for already done files and update the output set. """

        print("En _checkNewOutput")
        sys.stdout.flush()

        # Load previously done items (from text file)
        doneList = self._readDoneList()
        print("len(doneList)",len(doneList))
        sys.stdout.flush()

        # Check for newly done items
        particlesListId = self._readParticlesId()
        print("len(particlesListId)", len(particlesListId))
        sys.stdout.flush()

        newDone = [particlesId for particlesId in particlesListId
                           if particlesId not in doneList]
        #firstTime = len(doneList) == 0
        allDone = len(doneList) + len(newDone)

        #print("newDone, firstTime, allDone", newDone, firstTime, allDone)
        #sys.stdout.flush()

        # We have finished when there is not more inputs (stream closed)
        # and the number of processed particles is equal to the number of inputs
        self.finished = (self.isStreamClosed == Set.STREAM_CLOSED
                         and allDone == len(self.listOfParticles))
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if newDone:
            for particleId in newDone:
                self._writeDoneList(particleId)
        #elif not self.finished:
        #    print("Noooo")
        #    sys.stdout.flush()
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
        #    return

        outSet = self._loadOutputSet(SetOfClasses2D, 'classes2D.sqlite')

        if (exists(self._getPath('classes2D.sqlite'))):
            if(exists(self._getExtraPath('last_classes.xmd')) and exists(self._getExtraPath('last_images.xmd'))):
                self._updateOutputSet('outputParticles', outSet, streamMode)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

        if (exists(self._getPath('classes2D.sqlite'))):
            outSet.close()


    def _generateInputMd(self, reclassification):
        print("En _generateInputMd")
        sys.stdout.flush()

        doneList = self._readDoneList()
        metadataItem = md.MetaData(self.imgsExp)
        metadataInput = md.MetaData()

        #print("self.insertedListFlag",self.insertedListFlag)
        rows = iterRows(metadataItem)
        prevObjId = -1
        for row in rows:
            objId = row.getValue(md.MDL_ITEM_ID)
            if objId in self.insertedList:
                if objId not in doneList:
                    if prevObjId != -1:
                        if self.insertedListFlag[objId]!=self.insertedListFlag[prevObjId]:
                            break
                    prevObjId = objId
                    row.addToMd(metadataInput)


        metadataInput.write(self._getExtraPath('inputImagesExp.xmd'), MD_APPEND)
        print("metadataInput.len",getSize(self._getExtraPath('inputImagesExp.xmd')))

        fn = self._getExtraPath('particles.txt')
        rows = iterRows(metadataInput)
        for row in rows:
            with open(fn, 'a') as f:
                f.write('%d\n' % row.getValue(md.MDL_ITEM_ID))

        if reclassification:
            lenRefs = len(self.listNumImgs)
            randRef = randint(1, lenRefs)
            print("RECLASSIFICATION", randRef)
            outMd = self._getExtraPath('last_classes.xmd')
            block = 'class%06d' % (randRef)
            self._params = {'newMd': block + "_images@" + outMd,
                            'outMd': self._getExtraPath('inputImagesExp.xmd'),
                            }
            args = ('-i %(newMd)s -o %(outMd)s --set union_all %(outMd)s')
            self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)

            self._params = {'newMd': self._getExtraPath('last_images.xmd'),
                            'label': randRef,
                            }
            args = ('-i %(newMd)s -o %(newMd)s --query select "ref != %(label)d"')
            self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)

            mdAll = md.MetaData()

            fn = self._getExtraPath('last_classes.xmd')
            blocks = md.getBlocksInMetaDataFile(fn)
            count=1
            for block in blocks:
                if block.startswith('classes'):
                    mdClass = md.MetaData(block + "@" + fn)
                    rows = iterRows(mdClass)
                    for row in rows:
                        if mdClass.getValue(md.MDL_REF, row.getObjId()) == randRef:
                            row.setValue(md.MDL_CLASS_COUNT, 0L)
                        row.addToMd(mdAll)
            mdAll.write('classes@' + fn, MD_APPEND)
            for block in blocks:
                mdAll2 = md.MetaData()
                if block.startswith('class00'):
                    mdClass = md.MetaData(block + "@" + fn)
                    if not block.startswith('class%06d' % (randRef)):
                        rows = iterRows(mdClass)
                        for row in rows:
                            row.addToMd(mdAll2)
                    mdAll2.write('class%06d' % (count) + '_images@' + fn, MD_APPEND)
                    count+=1
            copy(self._getExtraPath('last_classes.xmd'), self._getExtraPath('prueba_classes.xmd'))
            copy(self._getExtraPath('last_images.xmd'), self._getExtraPath('prueba_images.xmd'))





    def generateMetadata(self, listNameImgs, listNumImgs, listRefImgs, level):

        print('generateMetadata', level)
        print('listNameImgs en generateMetadata', listNameImgs)
        print('listNumImgs en generateMetadata', listNumImgs)

        mdAll = md.MetaData()
        #mdImages = md.MetaData()
        count = 1

        for i in range(len(listNumImgs)):
            if listNumImgs[i] is not -1:
                name = listNameImgs[i]
                numRef = listRefImgs[i]
                print(name)

                if exists(self._getExtraPath(join('level%03d' % (level-1),'intermediate_classes.xmd'))):
                    fn = self._getExtraPath(join('level%03d' % (level-1),'intermediate_classes.xmd'))
                else:
                    fn = name[name.find('@')+1:-4]+'.xmd'
                blocks = md.getBlocksInMetaDataFile(fn)
                #if fn.find('final') != -1:
                #    fnImages = fn[:-17]+'final_images.xmd'
                #if fn.find('level_00') != -1:
                #    fnImages = fn[:-26]+'images.xmd'

                for block in blocks:
                    if block.startswith('classes'):
                        mdClass = md.MetaData(block + "@" + fn)
                        rows = iterRows(mdClass)
                        for row in rows:
                            if mdClass.getValue(md.MDL_REF, row.getObjId()) == numRef:

                                #args = ('-i %(newMd)s --query select "ref == %(numRef)d" -o %(outMd)s')
                                #self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)
                                #args2 = ('-i %(outMd)s --operate modify_values "ref = %(newNumRef)d" -o %(outMd)s')
                                #self.runJob("xmipp_metadata_utilities", args2 % self._params, numberOfMpi=1)
                                #print("getSize(outImages.xmd)", getSize(self._getExtraPath(join('level%03d' % level, 'outImages.xmd'))))
                                #mdOrig = md.MetaData(self._getExtraPath(join('level%03d' % level, 'outImages.xmd')))
                                #rowsO = iterRows(mdOrig)
                                #for rowO in rowsO:
                                #    rowO.addToMd(mdImages)
                                #mdImages.write(self._getExtraPath(join('level%03d' % level, 'intermediate_images.xmd')),  MD_APPEND)

                                row.setValue(md.MDL_REF, count)
                                row.addToMd(mdAll)
                                count = count + 1
                                break

        outSet = self._getExtraPath(join('level%03d' % level,'level%03d' % level + '_classes.xmd'))
        blocks = md.getBlocksInMetaDataFile(outSet)
        #i=1
        for block in blocks:
            if block.startswith('classes'):
                mdClass = md.MetaData(block + "@" + outSet)
                rows = iterRows(mdClass)
                for row in rows:

                    #args = ('-i %(newMd)s --query select "ref == %(numRef)d" -o %(outMd)s')
                    #self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)
                    #args2 = ('-i %(outMd)s --operate modify_values "ref = %(newNumRef)d" -o %(outMd)s')
                    #self.runJob("xmipp_metadata_utilities", args2 % self._params, numberOfMpi=1)
                    #print("getSize(outImages.xmd)", getSize(self._getExtraPath(join('level%03d' % level, 'outImages.xmd'))))
                    #mdOrig = md.MetaData(self._getExtraPath(join('level%03d' % level, 'outImages.xmd')))
                    #rowsO = iterRows(mdOrig)
                    #for rowO in rowsO:
                    #    rowO.addToMd(mdImages)
                    #mdImages.write(self._getExtraPath(join('level%03d' % level, 'intermediate_images.xmd')), MD_APPEND)

                    row.setValue(md.MDL_REF, count)
                    row.addToMd(mdAll)
                    count = count + 1
                    #i+=1

        mdAll.write('classes@' + self._getExtraPath(join('level%03d' % level,'intermediate_classes.xmd')), MD_APPEND)

        count = 1
        for i in range(len(listNumImgs)):
            if listNumImgs[i] is not -1:
                name = listNameImgs[i]
                numRef = listRefImgs[i]

                if exists(self._getExtraPath(join('level%03d' % (level-1),'intermediate_classes.xmd'))):
                    fn = self._getExtraPath(join('level%03d' % (level-1),'intermediate_classes.xmd'))
                else:
                    fn = name[name.find('@')+1:-4]+'.xmd'
                mdAll2 = md.MetaData()
                blocks = md.getBlocksInMetaDataFile(fn)

                for block in blocks:
                    if block.startswith('class%06d' % (numRef)):

                        self._params = {'newMd': block + "@" + fn,
                                        'outMd': self._getExtraPath(join('level%03d' % level, 'outImages.xmd')),
                                        'finalMd': self._getExtraPath(join('level%03d' % level, 'intermediate_images.xmd')),
                                        'newNumRef': count}
                        if count == 1:
                            args = ('-i %(newMd)s -o %(finalMd)s --operate modify_values "ref=%(newNumRef)d"')
                            self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)
                        else:
                            args = ('-i %(newMd)s -o %(outMd)s --operate modify_values "ref=%(newNumRef)d"')
                            self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)
                            args = ('-i %(outMd)s -o %(finalMd)s --set union_all %(finalMd)s')
                            self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)


                        mdClass = md.MetaData(block + "@" + fn)
                        rows = iterRows(mdClass)
                        for row in rows:
                            row.setValue(md.MDL_REF, count)
                            row.addToMd(mdAll2)
                        mdAll2.write('class%06d' % (count) + '_images@' + self._getExtraPath(join('level%03d' % level,'intermediate_classes.xmd')), MD_APPEND)
                        count = count + 1
                        break

        outSet = self._getExtraPath(join('level%03d' % level,'level%03d' % level + '_classes.xmd'))
        blocks = md.getBlocksInMetaDataFile(outSet)
        for block in blocks:
            mdAll2 = md.MetaData()
            if block.startswith('class0'):

                self._params = {'newMd': block + "@" + outSet,
                                'outMd': self._getExtraPath(join('level%03d' % level, 'outImages.xmd')),
                                'finalMd': self._getExtraPath(join('level%03d' % level, 'intermediate_images.xmd')),
                                'newNumRef': count}
                if count == 1:
                    args = ('-i %(newMd)s -o %(finalMd)s --operate modify_values "ref=%(newNumRef)d"')
                    self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)
                else:
                    args = ('-i %(newMd)s -o %(outMd)s --operate modify_values "ref=%(newNumRef)d"')
                    self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)
                    args = ('-i %(outMd)s -o %(finalMd)s --set union_all %(finalMd)s')
                    self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)


                mdClass = md.MetaData(block + "@" + outSet)
                rows = iterRows(mdClass)
                for row in rows:
                    row.setValue(md.MDL_REF, count)
                    row.addToMd(mdAll2)
                mdAll2.write('class%06d' % (count) + '_images@' + self._getExtraPath(join('level%03d' % level,'intermediate_classes.xmd')), MD_APPEND)
                count = count + 1

        self.iterReturnClass = 0

    def _splitStep(self, expImgMd, iterReturnSplit, flag_attraction):

        if getSize(expImgMd) == 0:
            print("NOoooooooooooooooooooooooo")
            return

        level = self.level
        iterReturnSplit = self.splitStep(expImgMd, level, iterReturnSplit, flag_attraction)
        #############################################################
        if not self.useAttraction:
            self.attractionSplitStep(level, iterReturnSplit)
        ##############################################################
        self.generateMetadata(self.listNameImgs, self.listNumImgs, self.listRefImgs, level)
        copy(self._getExtraPath(join('level%03d' % level, 'intermediate_classes.xmd')),
             self._getExtraPath('last_classes.xmd'))
        copy(self._getExtraPath(join('level%03d' % level, 'intermediate_images.xmd')),
             self._getExtraPath('last_images.xmd'))


    def _classifyStep(self, expImgMd, iterReturnClass):

        if getSize(expImgMd) == 0:
            print("Noooooooooooooooooooooooooo")
            return

        level = self.level
        refImgMd = self._getExtraPath('last_classes.xmd')
        iterReturnClass = self.classifyWholeSetStep(refImgMd, expImgMd, level, iterReturnClass)
        ##############################################################################
        if not self.useAttraction:
            change = True
            while change:
                change = self.attractionGeneralStep(level)
                if change:
                    level = level + 1
                    self._splitStep(expImgMd, level, 0, False)
                    self._classifyStep(expImgMd, level, iterReturnClass)
        ##############################################################################
        #copy(self._getExtraPath(join('level%03d' % level, 'general_level%03d' % level + '_classes.xmd')),
        #     self._getExtraPath('last_classes.xmd'))
        #copy(self._getExtraPath(join('level%03d' % level, 'general_images_level%03d' % level + '.xmd')),
        #     self._getExtraPath('last_images.xmd'))
        self.generateOutputClasses(level)


    def averageClasses(self, finalMetadata, lastMetadata, newMetadata, flag_iter):

        listLastClasses = []
        listFnLastClasses = []
        listNewClasses = []
        listFnNewClasses = []
        listRefLast = []

        metadataItemLast = md.MetaData(lastMetadata)
        metadataItemNew = md.MetaData(newMetadata)
        if flag_iter:
            metadataFinal = md.MetaData(finalMetadata)
            nameFinal = metadataFinal.getValue(md.MDL_IMAGE, 1)
            nameFinal = nameFinal[7:-3]+'stk'
        else:
            nameFinal = finalMetadata
        print("nameFinal",nameFinal)

        for item in metadataItemLast:
            numImgs = metadataItemLast.getValue(md.MDL_CLASS_COUNT, item)
            listLastClasses.append(numImgs)
            nameRef = metadataItemLast.getValue(md.MDL_IMAGE, item)
            listFnLastClasses.append(nameRef)
            listRefLast.append(metadataItemLast.getValue(md.MDL_REF, item))

        print("listRefLast", listRefLast)
        i=0
        for item in metadataItemNew:
            labelRefNew = metadataItemNew.getValue(md.MDL_REF, item)
            nameRef = metadataItemNew.getValue(md.MDL_IMAGE, item)
            print("labelRefNew",labelRefNew)
            if labelRefNew>listRefLast[i]:
                print("NOOOOOO")
                #listNewClasses.append(0)
                #listFnNewClasses.append('None')
                #i+=1
            numImgs = metadataItemNew.getValue(md.MDL_CLASS_COUNT, item)
            listNewClasses.append(numImgs)
            listFnNewClasses.append(nameRef)
            i+=1


        print("listFnLastClasses", listFnLastClasses)
        print("listFnNewClasses", listFnNewClasses)
        print("listLastClasses", listLastClasses)
        print("listNewClasses", listNewClasses)
        total = []
        for i in range(len(listLastClasses)):
            listToMultiply = []
            total.append(listLastClasses[i] + listNewClasses[i])
            listToMultiply.append(float(listLastClasses[i]) / float(total[i]))
            listToMultiply.append(float(listNewClasses[i]) / float(total[i]))
            print("total", total)
            print("listToMultiply", listToMultiply)

            if listNewClasses[i]==0:
                im1 = xmipp.Image(listFnLastClasses[i])
                im1.write('%06d@' % (i + 1) + nameFinal)
                continue
            if listLastClasses[i]==0:
                im2 = xmipp.Image(listFnNewClasses[i])
                im2.write('%06d@' % (i + 1) + nameFinal)
                continue

            im1 = xmipp.Image(listFnLastClasses[i])
            im2 = xmipp.Image(listFnNewClasses[i])
            im1.inplaceMultiply(listToMultiply[0])
            im2.inplaceMultiply(listToMultiply[1])
            im1.align(im2)

            im1.convert2DataType(xmipp.DT_DOUBLE)
            im2.convert2DataType(xmipp.DT_DOUBLE)

            im1.inplaceAdd(im2)  # Aligned

            im1.write('%06d@' % (i + 1) + nameFinal)

            # self._params = {'imgInOne': listFnLastClasses[i],
            #                'imgOutOne': self._getExtraPath('classOne.stk'),
            #                'valueOne': listToMultiply[0]}
            # args = ('-i %(imgInOne)s --mult %(valueOne)f -o %(imgOutOne)s')
            # self.runJob("xmipp_image_operate", args % self._params, numberOfMpi=1)

            # self._params = {'imgInTwo': listFnNewClasses[i],
            #                'imgOutTwo': self._getExtraPath('classTwo.stk'),
            #                'valueTwo': listToMultiply[1]}
            # args = ('-i %(imgInTwo)s --mult %(valueTwo)f -o %(imgOutTwo)s')
            # self.runJob("xmipp_image_operate", args % self._params, numberOfMpi=1)

            # self._params = {'imgInOne': self._getExtraPath('alignedOne.stk'),
            #                'imgInTwo': self._getExtraPath('alignedTwo.stk'),
            #                'imgOut': '%06d@'%(i+1) + finalMetadata }
            # args = ('-i %(imgInOne)s --plus %(imgInTwo)s -o %(imgOut)s')
            # self.runJob("xmipp_image_operate", args % self._params, numberOfMpi=1)

        if flag_iter:
            copy(finalMetadata, self._getExtraPath('prueba_average%06d.stk')%self.level)
        return total


    def generateOutputClasses(self, level):

        print("En generateOutputClasses")
        if level==1:
            copy(self._getExtraPath(join('level%03d' % level, 'general_level%03d' % level + '_classes.xmd')),
                 self._getExtraPath('last_classes.xmd'))
            copy(self._getExtraPath(join('level%03d' % level, 'general_images_level%03d' % level + '.xmd')),
                 self._getExtraPath('last_images.xmd'))
            return

        finalMetadata = self._getExtraPath('final_classes.stk')
        lastMetadata = self._getExtraPath('last_classes.xmd')
        newMetadata = self._getExtraPath(join('level%03d' % level, 'general_level%03d' % level + '_classes.xmd'))

        total = self.averageClasses(finalMetadata, lastMetadata, newMetadata, False)

        copy(self._getExtraPath('final_classes.stk'), self._getExtraPath('last_classes.stk'))

        mdAll = md.MetaData()
        for i in range(len(total)):
            row = md.Row()
            row.setValue(md.MDL_REF, i+1)
            row.setValue(md.MDL_IMAGE, '%06d@'%(i+1) + finalMetadata)
            row.setValue(md.MDL_CLASS_COUNT, total[i])
            row.addToMd(mdAll)
        mdAll.write('classes@'+finalMetadata[:-3] + 'xmd', MD_APPEND)

        for i in range(len(total)):
            self._params = {'lastMd': 'class%06d_images@' % (i + 1) + lastMetadata,
                            'newMd': 'class%06d_images@' % (i + 1) + newMetadata,
                            'outMd': 'class%06d_images@' % (i + 1) + finalMetadata[:-3] + 'xmd'}
            args = ('-i %(lastMd)s --set union_all %(newMd)s -o %(outMd)s --mode append')
            self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)

        copy(self._getExtraPath('final_classes.xmd'), self._getExtraPath('last_classes.xmd'))

        self._params = {'lastMd': self._getExtraPath('last_images.xmd'),
                        'newMd': self._getExtraPath(join('level%03d' % level, 'general_images_level%03d' % level + '.xmd')),
                        'outMd': self._getExtraPath('final_images.xmd')}
        args = ('-i %(lastMd)s --set union %(newMd)s -o %(outMd)s')
        self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)

        copy(self._getExtraPath('final_images.xmd'), self._getExtraPath('last_images.xmd'))



    def splitStep(self, expImgMd, level, iterReturnSplit, flag_attraction):

        print('splitStep')

        i=0
        refImgMd=None
        while i <self.numberOfSplitIterations:
            self.iterationStep(refImgMd, expImgMd, i, True, flag_attraction, level)
            refImgMd = self._getExtraPath(join('level%03d' % level,'level%03d' % level + '_classes.xmd'))
            i+=1

            # AJ comprobar si hay dos clases o no porque si no las hay, TENEMOS UN PROBLEMA
            length = getSize(self._getExtraPath(join('level%03d' % level,'level%03d' % level + '_classes.xmd')))
            if length == 1:
                print('PROBLEEEEEEM')
                i = 0

            if not self.useAttraction:
                if self.fastCheckingAtt(level, True):
                    iterReturnSplit=i
                    return iterReturnSplit

        return iterReturnSplit


    def classifyWholeSetStep(self, refImgMd, expImgMd, level, iterReturnClass):

        print('classifyWholeSetStep')

        i=iterReturnClass
        if i == self.numberOfClassifyIterations:
            i -= 1
        while i <self.numberOfClassifyIterations:
            print("ITER",i)
            self.iterationStep(refImgMd, expImgMd, i, False, False, level)

            refImgMd = self._getExtraPath(join('level%03d' % level,'general_level%03d' % level + '_classes.xmd'))

            if self.checkContinueClassification(level, i):
                return

            #AJ NEW
            if i+1 < self.numberOfClassifyIterations:
                finalMetadata = self._getExtraPath(join('level%03d' % level,'general_level%03d' % level + '_classes.xmd'))
                #lastMetadata = self._getExtraPath('last_classes.xmd')
                if exists(self._getExtraPath('final_classes.xmd')):
                    lastMetadata = self._getExtraPath('final_classes.xmd')
                else:
                    if level<2:
                        lastMetadata = self._getExtraPath('last_classes.xmd')
                    else:
                        lastMetadata = self._getExtraPath(join('level%03d' % (level-1), 'general_level%03d' % (level-1) + '_classes.xmd'))
                newMetadata = self._getExtraPath(join('level%03d' % level, 'general_level%03d' % level + '_classes.xmd'))
                print("lastMetadata", i, lastMetadata)
                print("finalMetadata", finalMetadata)
                print("newMetadata", newMetadata)
                self.averageClasses(finalMetadata, lastMetadata, newMetadata, True)
            #FIN AJ

            i += 1

            # AJ check attraction
            if not self.useAttraction:
                if self.fastCheckingAtt(level, False):
                    iterReturnClass=i
                    return iterReturnClass
        return iterReturnClass


    # AJ check attraction
    def fastCheckingAtt (self, level, flag_split):

        if flag_split:
            metadata = md.MetaData(self._getExtraPath(join('level%03d' % level, 'level%03d' % level + '_classes.xmd')))
            total = getSize(self._getExtraPath(join('level%03d' % level, 'images_level%03d' % level + '.xmd')))
            print("total", total)
            th = (self.p * total / 2)
            for item in metadata:
                numImgs = metadata.getValue(md.MDL_CLASS_COUNT, item)
                print("numImgs", numImgs)
                print("th", th)
                if numImgs < th:
                    return True
            return False
        else:
            listAuxNum = []
            metadata = md.MetaData(self._getExtraPath(join('level%03d' % level, 'general_level%03d' % level + '_classes.xmd')))
            for item in metadata:
                numImgs = metadata.getValue(md.MDL_CLASS_COUNT, item)
                listAuxNum.append(numImgs)
                print("numImgs",numImgs)
            total = sum(listAuxNum)
            th = (self.p * total / len(listAuxNum))
            print("th", th)
            aux = [i for i in listAuxNum if i<th]
            if len(aux)>0:
                return True
            else:
                return False

    def checkContinueClassification(self, level, iter):

        print('checkContinueClassification', level, iter)
        diff=0
        i=0
        metadata = md.MetaData(self._getExtraPath(join('level%03d' % level, 'general_images_level%03d' % level + '.xmd')))

        for item in metadata:
            refImg = metadata.getValue(md.MDL_REF, item)
            nameImg = metadata.getValue(md.MDL_IMAGE, item)
            if iter==0:
                self.listContinueClass.append(nameImg)
                self.listContinueClass.append(refImg)
            else:
                if nameImg in self.listContinueClass:
                    idx = self.listContinueClass.index(nameImg) + 1
                    if refImg!=self.listContinueClass[idx]:
                        diff+=1
                    self.listContinueClass[idx]=refImg
                else:
                    diff += 1
                    self.listContinueClass.append(nameImg)
                    self.listContinueClass.append(refImg)
            i+=1
        num=(diff*100/i)
        print('checkContinueClassification',num,diff,i)
        if num<self.percentStopClassify and iter>0:
            return True
        else:
            return False


    def iterationStep (self, refSet, imgsExp, iter, flag_split, flag_attraction, level):

        print('iterationStep')

        if not exists(join(self._getExtraPath(), 'level%03d' % level)):
            mkdir(join(self._getExtraPath(), 'level%03d' % level))

        if iter==0 and flag_split==True:

            # First step: divide the metadata input file to generate
            # a couple of references
            if level==0:
                if not flag_attraction:
                    outDirName = imgsExp[0:imgsExp.find('extra')+6] + 'level%03d' % level + imgsExp[imgsExp.find('extra')+5:-4]
                else:
                    outDirName = imgsExp[0:imgsExp.find('extra') + 6] + 'level%03d' % level + imgsExp[imgsExp.find('level%03d' % level) + 8:-4]
            else:
                if not flag_attraction:
                    outDirName = imgsExp[0:imgsExp.find('extra') + 6] + 'level%03d' % level + imgsExp[imgsExp.find('level%03d' % (level-1)) + 8:-4]
                else:
                    outDirName = imgsExp[0:imgsExp.find('extra') + 6] + 'level%03d' % level + imgsExp[imgsExp.find('level%03d' % level) + 8:-4]
            self._params = {'imgsExp': imgsExp,
                            'outDir': outDirName}
            args = ('-i %(imgsExp)s -n 2 --oroot %(outDir)s')
            self.runJob("xmipp_metadata_split", args % self._params, numberOfMpi=1)

            # Second step: calculate the means of the previous metadata
            expSet1 = outDirName + '000001.xmd'
            avg1 = outDirName + '_000001'
            expSet2 = outDirName + '000002.xmd'
            avg2 = outDirName + '_000002'
            self._params = {'imgsSet': expSet1,
                            'outputAvg': avg1}
            args = ('-i %(imgsSet)s --save_image_stats %(outputAvg)s -v 0')
            self.runJob("xmipp_image_statistics", args % self._params, numberOfMpi=1)

            self._params = {'imgsSet': expSet2,
                            'outputAvg': avg2}
            args = ('-i %(imgsSet)s --save_image_stats %(outputAvg)s -v 0')
            self.runJob("xmipp_image_statistics", args % self._params, numberOfMpi=1)

            # Third step: generate a single metadata with the two previous averages
            refSet = self._getExtraPath(join('level%03d' % level,'refSet.xmd'))
            self._params = {'avg1': avg1 + 'average.xmp',
                            'avg2': avg2 + 'average.xmp',
                            'outputMd': refSet}
            args = ('-i %(avg1)s --set union %(avg2)s -o %(outputMd)s')
            self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)

        # Fourth step: calling program xmipp_cuda_correlation
        if flag_split:
            filename = 'level%03d' % level+'_classes.xmd'
            self._params = {'imgsRef': refSet,
                            'imgsExp': imgsExp,
                            'outputFile': 'images_level%03d' % level+'.xmd',
                            'tmpDir': join(self._getExtraPath(),'level%03d' % level),
                            'keepBest': self.keepBest.get(),
                            'maxshift': self.maximumShift.get(),
                            'outputClassesFile': filename,
                            }
        else:
            filename = 'general_level%03d' % level + '_classes.xmd'
            self._params = {'imgsRef': refSet,
                            'imgsExp': imgsExp,
                            'outputFile': 'general_images_level%03d' % level + '.xmd',
                            'tmpDir': join(self._getExtraPath(),'level%03d' % level),
                            'keepBest': self.keepBest.get(),
                            'maxshift': self.maximumShift.get(),
                            'outputClassesFile': filename,
                            }
        if not flag_split:
            args = '-i_ref %(imgsRef)s -i_exp %(imgsExp)s -o %(outputFile)s '\
                   '--odir %(tmpDir)s --keep_best %(keepBest)d '\
                   '--maxShift %(maxshift)d --classify %(outputClassesFile)s '\
                   '--simplifiedMd'
            self.runJob("xmipp_cuda_correlation", args % self._params, numberOfMpi=1)
        else:
            Nrefs = getSize(refSet)
            args = '-i %(imgsExp)s --ref0 %(imgsRef)s --nref %(Nrefs)d --iter 1 --distance correlation '\
                   '--classicalMultiref --maxShift %(maxshift)d --odir %(cl2dDir)s'
            self._params['Nrefs']=Nrefs
            self._params['cl2dDir'] = self._getExtraPath(join('level%03d' % level))
            self.runJob("xmipp_classify_CL2D", args % self._params)
            copy(self._getExtraPath(join('level%03d' % level,"level_00","class_classes.xmd")),
                 self._getExtraPath(join('level%03d' % level,'level%03d' % level+'_classes.xmd')))
            copy(self._getExtraPath(join('level%03d' % level,"images.xmd")),
                 self._getExtraPath(join('level%03d' % level,'images_level%03d' % level + '.xmd')))

        #AJ only CL2D
        #if flag_split:
        #    copy(self._getExtraPath(join('level%03d' % level, "level_00", "class_classes.xmd")),
        #         self._getExtraPath(join('level%03d' % level, 'level%03d' % level + '_classes.xmd')))
        #    copy(self._getExtraPath(join('level%03d' % level, "images.xmd")),
        #         self._getExtraPath(join('level%03d' % level, 'images_level%03d' % level + '.xmd')))
        #else:
        #    copy(self._getExtraPath(join('level%03d' % level, "level_00", "class_classes.xmd")),
        #         self._getExtraPath(join('level%03d' % level, 'general_level%03d' % level + '_classes.xmd')))
        #    copy(self._getExtraPath(join('level%03d' % level, "images.xmd")),
        #         self._getExtraPath(join('level%03d' % level, 'general_images_level%03d' % level + '.xmd')))




    def attractionSplitStep(self, level, iterReturnSplit):

        change, labelMaxClass, labelMinClass, mdToReduce, self.listNumImgs, self.listNameImgs = self.checkAttraction(level, True)

        while change:

            self.depthSplit += 1

            print('CHANGEEEEEEEEEEEEEEEE')

            self._params = {'input': mdToReduce,
                            'numRef': labelMaxClass,
                            'outputMd': mdToReduce[0:-4]+'_NoAtt.xmd'}
            args = ('-i %(input)s --query select "ref==%(numRef)d" -o %(outputMd)s')
            self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)


            self.imgsExp = self._getExtraPath(join('level%03d' % level,'images_level%03d' % level + '_NoAtt.xmd'))

            i = 0
            while i < self.numberOfSplitIterations:
                self.iterationStep(self.refSet, self.imgsExp, i, True, True, level)
                self.refSet = self._getExtraPath(join('level%03d' % level,'level%03d' % level + '_classes.xmd'))
                i+=1

                # AJ comprobar si hay dos clases o no porque si no las hay, TENEMOS UN PROBLEMA
                length = getSize(self._getExtraPath(join('level%03d' % level,'level%03d' % level + '_classes.xmd')))
                if length == 1:
                    print('PROBLEEEEEEM')
                    i = 0

                if self.fastCheckingAtt(level, True):
                    break

            self.attractionSplitStep(level)

            if self.depthSplit>1:
                self.depthSplit-=1
                return
            if self.depthSplit==1:
                self.depthSplit=0

            if (level - 1) >= 0:
                self.imgsExp = self._getExtraPath(join('level%03d' % (level-1),'images_level%03d' % (level-1) + '_major.xmd'))
            else:
                self.imgsExp = self._getExtraPath('imagesExp.xmd')

            i = iterReturnSplit
            if i==self.numberOfSplitIterations:
                i-=1
            while i < self.numberOfSplitIterations:
                self.iterationStep(self.refSet, self.imgsExp, i, True, False, level)
                self.refSet = self._getExtraPath(join('level%03d' % level,'level%03d' % level + '_classes.xmd'))
                i+=1

                # AJ comprobar si hay dos clases o no porque si no las hay, TENEMOS UN PROBLEMA
                length = getSize(self._getExtraPath(join('level%03d' % level,'level%03d' % level + '_classes.xmd')))
                if length == 1:
                    print('PROBLEEEEEM')
                    i = 0

                if self.fastCheckingAtt(level, True):
                    break

            change, labelMaxClass, __, mdToReduce, __, __ = self.checkAttraction(level, True)



    def attractionGeneralStep(self, level):

        change, labelMaxClass, labelMinClass, mdToReduce, self.listNumImgs, self.listNameImgs = self.checkAttraction(level, False)

        if change:
            self.depth+=1

        if change: #esto era un while
            print('CHANGEEEEEEEEEEEEEEEE', level)
            print('listNameImgs', self.listNameImgs)
            print('listNumImgs', self.listNumImgs)
            print('mdToReduce', mdToReduce)
            print('labelMaxClass', labelMaxClass)

            self._params = {'input': mdToReduce,
                            'numRef': labelMaxClass,
                            #'outputMd': mdToReduce[0:-4]+'_NoAtt.xmd'}
                            'outputMd': mdToReduce[0:-25] + 'images_level%03d' % level + '_major.xmd'}
            args = ('-i %(input)s --query select "ref==%(numRef)d" -o %(outputMd)s')
            self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)

        return change


    def checkAttraction(self, level, flag_split):

        if flag_split:
            mdToCheck = self._getExtraPath(join('level%03d' % level,'level%03d' % level+'_classes.xmd'))
            mdToReduce = self._getExtraPath(join('level%03d' % level,'images_level%03d' % level+'.xmd'))
        else:
            mdToCheck = self._getExtraPath(join('level%03d' % level,'general_level%03d' % level + '_classes.xmd'))
            mdToReduce = self._getExtraPath(join('level%03d' % level,'general_images_level%03d' % level + '.xmd'))

        listAuxNum=[]
        listAuxName = []
        listAuxRef = []
        metadataItem = md.MetaData(mdToCheck)
        for item in metadataItem:
            numImgs = metadataItem.getValue(md.MDL_CLASS_COUNT, item)
            name = metadataItem.getValue(md.MDL_IMAGE, item)
            ref = metadataItem.getValue(md.MDL_REF, item)
            listAuxNum.append(numImgs)
            listAuxName.append(name)
            listAuxRef.append(ref)
            print(listAuxNum)
        total = sum(listAuxNum)

        th = (self.p*total/len(listAuxNum))
        labelMinClass=[]
        print('TH', self.p * total / len(listAuxNum))
        for i in range(len(listAuxNum)):
            if listAuxNum[i]<th:
                labelMinClass.append(listAuxRef[i])
                listAuxNum[i] = -1

        labelMaxClass = listAuxRef[listAuxNum.index(max(listAuxNum))]
        print("labelMaxClass",labelMaxClass)
        listAuxNum[labelMaxClass-1] = -1

        if len(labelMinClass)>0:
            change = True
        else:
            change = False

        return change, labelMaxClass, labelMinClass, mdToReduce, listAuxNum, listAuxName


    def checkSplit(self, expImgMd):

        if getSize(expImgMd) == 0:
            print("NOoooooooooooooooooooooooooooooooooo")
            return

        print("In checkSplit")

        outSet = self._getExtraPath('last_classes.xmd')

        level = self.level
        self.listNameImgs = []
        self.listNumImgs = []
        self.listRefImgs = []

        metadataItem = md.MetaData(outSet)
        for item in metadataItem:
            nameImg = metadataItem.getValue(md.MDL_IMAGE, item)
            self.listNameImgs.append(nameImg)
            numImgs = metadataItem.getValue(md.MDL_CLASS_COUNT, item)
            self.listNumImgs.append(numImgs)
            refImg = metadataItem.getValue(md.MDL_REF, item)
            self.listRefImgs.append(refImg)
        print(self.listNameImgs)
        print(self.listNumImgs)
        print(self.listRefImgs)

        total = sum(self.listNumImgs)

        thPercent = (1.2 * total / len(self.listNumImgs))
        thGlobal = 500
        i=0

        while i<len(self.listNumImgs):

            print("listAuxNum", self.listNumImgs[i])
            print("total", total)
            print("TH", thPercent)

            #if self.listNumImgs[i]<thPercent and self.listNumImgs[i]<thGlobal:
            if self.listNumImgs[i]<1.75*thGlobal:
                i+=1
                continue

            maxValue = self.listNumImgs[i]
            maxPos = self.listNumImgs.index(maxValue)

            self.listNumImgs[maxPos] = -1
            bestRef = self.listRefImgs[maxPos]

            outputMd = self._getExtraPath(join('level%03d' % level,'general_images_level%03d' % level + '_major.xmd'))
            self._params = {'input': 'class%06d_images' % (bestRef) + '@' + outSet,
                            'outputMd': outputMd
                            }
            args = ('-i %(input)s -o %(outputMd)s')
            self.runJob("xmipp_metadata_utilities", args % self._params, numberOfMpi=1)

            #Split with XXXX_major.xmd
            #Generate metadata with the output of the previous split and the data in general_XXX_classes.xmd (removing the major class)
            self.levelUp(outputMd)
            if not exists(join(self._getExtraPath(),'level%03d' % self.level)):
                mkdir(join(self._getExtraPath(),'level%03d' % self.level))
            newOutputMd = self._getExtraPath(join('level%03d' % self.level, 'general_images_level%03d' % self.level + '_major.xmd'))
            copy(outputMd, newOutputMd)
            self._splitStep(newOutputMd, 0, True)

            i=0
            outSet = self._getExtraPath('last_classes.xmd')

            level = self.level
            self.listNameImgs = []
            self.listNumImgs = []
            self.listRefImgs = []

            metadataItem = md.MetaData(outSet)
            for item in metadataItem:
                nameImg = metadataItem.getValue(md.MDL_IMAGE, item)
                self.listNameImgs.append(nameImg)
                numImgs = metadataItem.getValue(md.MDL_CLASS_COUNT, item)
                self.listNumImgs.append(numImgs)
                refImg = metadataItem.getValue(md.MDL_REF, item)
                self.listRefImgs.append(refImg)
            print(self.listNameImgs)
            print(self.listNumImgs)
            print(self.listRefImgs)

            copy(outSet, self._getExtraPath('final_classes.xmd'))


    def checkOutput(self, level):

        print('checkOutput')

        listAuxString = []
        listAuxNum = []
        listAuxRefs = []

        outSet = self._getExtraPath(join('level%03d' % level,'intermediate_classes.xmd'))
        metadataItem = md.MetaData(outSet)
        for item in metadataItem:
            nameImg = metadataItem.getValue(md.MDL_IMAGE, item)
            listAuxString.append(nameImg)
            numImgs = metadataItem.getValue(md.MDL_CLASS_COUNT, item)
            listAuxNum.append(numImgs)
            refImg = metadataItem.getValue(md.MDL_REF, item)
            listAuxRefs.append(refImg)
        print(listAuxString)
        print(listAuxNum)
        print(listAuxRefs)

        maxValue = max(listAuxNum)
        maxPos = listAuxNum.index(maxValue)

        listAuxNum[maxPos] = -1
        bestRef = listAuxRefs[maxPos]

        mdAll = md.MetaData()
        outSet = self._getExtraPath(join('level%03d' % level,'intermediate_classes.xmd'))
        blocks = md.getBlocksInMetaDataFile(outSet)
        for block in blocks:
            if block.startswith('class%06d' % (bestRef)):
                mdClass = md.MetaData(block + "@" + outSet)
                rows = iterRows(mdClass)
                for row in rows:
                    row.addToMd(mdAll)
        mdAll.write(self._getExtraPath(join('level%03d' % level,'images_level%03d' % level+'_major.xmd')), MD_APPEND)

        return listAuxNum, listAuxString


    def cleaningPath(self, level):
        if exists(self._getExtraPath(join('level%03d' % level,'refSet.xmd'))):
            remove(self._getExtraPath(join('level%03d' % level,'refSet.xmd')))
        if exists(self._getExtraPath(join('level%03d' % level,"images.xmd"))):
            remove(self._getExtraPath(join('level%03d' % level,"images.xmd")))
            level_1=level-1
        if level>0 and exists(self._getExtraPath(join('level%03d' % level_1,"images_level%03d_major.xmd" %level_1))):
            remove(self._getExtraPath(join('level%03d' % level_1,"images_level%03d_major.xmd" %level_1)))
        for f in listdir(join(join(self._getExtraPath(),'level%03d'%level))):
            if not f.find('average')==-1 or not f.find('stddev')==-1 or not f.find('NoAtt')==-1:
                remove(join(join(self._getExtraPath(),'level%03d'%level,f)))
            if level>0 and not f.find('major0') == -1:
                remove(join(join(self._getExtraPath(), 'level%03d' % level, f)))
            if level==0 and not f.find('imagesExp0') == -1:
                remove(join(join(self._getExtraPath(), 'level%03d' % level, f)))


    def createOutputStep(self):
        pass


    # --------------------------- UTILS functions -------------------------------

    def convertSet(self, imgsFileName, setOfImg):
        if setOfImg is None:
            return
        writeSetOfParticles(setOfImg, imgsFileName, alignType=ALIGN_NONE)


    def _loadInputList(self):
        """ Load the input set of ctfs and create a list. """
        particlesSet = self._loadInputParticleSet()
        self.isStreamClosed = particlesSet.getStreamState()
        self.listOfParticles = [m.getObjId() for m in particlesSet]
        particlesSet.close()
        self.debug("Closed db.")

    def _loadInputParticleSet(self):
        particlesFile = self.inputParticles.get().getFileName()
        self.debug("Loading input db: %s" % particlesFile)
        particlesSet = SetOfParticles(filename=particlesFile)
        particlesSet.loadAllProperties()
        return particlesSet

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == 'createOutputStep':
                return s
        return None

    def _readDoneList(self):
        """ Read from a text file the id's of the items that have been done. """
        DoneFile = self._getExtraPath('DONE.TXT')
        DoneList = []
        # Check what items have been previously done
        if exists(DoneFile):
            with open(DoneFile) as f:
                DoneList += [int(line.strip()) for line in f]
        return DoneList

    def _writeDoneList(self, particleId):
        """ Write to a text file the items that have been done. """
        AcceptedFile = self._getExtraPath('DONE.TXT')
        with open(AcceptedFile, 'a') as f:
            f.write('%d\n' % particleId)

    def _readParticlesId(self):

        fn = self._getExtraPath('particles.txt')
        particlesList = []
        # Check what items have been previously done
        if exists(fn):
            with open(fn) as f:
                particlesList += [int(line.strip().split()[0]) for line in f]
        return particlesList

    def _loadOutputSet(self, SetClass, baseName):
        setFile = self._getPath(baseName)
        if exists(setFile):
            outputSet = SetClass(filename=setFile)
            #outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        inputs = self.inputParticles.get()
        outputSet.copyInfo(inputs)
        return outputSet


    def _updateOutputSet(self, outputName, outputSet, state=Set.STREAM_OPEN):
        print("In _updateOutputSet")
        sys.stdout.flush()

        outputSet.setStreamState(state)
        if self.hasAttribute(outputName):
            print("In _updateOutputSet IF")
            sys.stdout.flush()

            outputSet.write()  # Write to commit changes
            outputAttr = getattr(self, outputName)
            # Copy the properties to the object contained in the protocol
            outputAttr.copy(outputSet, copyId=False)
            # Persist changes
            self._store(outputAttr)
        else:
            print("In _updateOutputSet ELSE")
            sys.stdout.flush()
            inputs = self.inputParticles.get()
            # Here the defineOutputs function will call the write() method
            outputSet = self._createSetOfClasses2D(inputs)
            #print("outputSet", outputSet.getImages())
            self._fillClassesFromLevel(outputSet)

            self._defineOutputs(**{'outputClasses': outputSet})
            self._defineSourceRelation(self.inputParticles, outputSet)
            self._store(outputSet)



    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(md.MDL_REF))
        item.setTransform(rowToAlignment(row, ALIGN_2D))


    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn, _ = self._classesInfo[classId]
            item.setAlignment2D()
            rep = item.getRepresentative()
            rep.setLocation(index, fn)
            rep.setSamplingRate(self.inputParticles.get().getSamplingRate())


    def _loadClassesInfo(self, filename):
        """ Read some information about the produced 2D classes
        from the metadata file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id

        mdClasses = md.MetaData(filename)
        for classNumber, row in enumerate(md.iterRows(mdClasses)):
            index, fn = xmippToLocation(row.getValue(md.MDL_IMAGE))
            self._classesInfo[classNumber + 1] = (index, fn, row.clone())



    def _fillClassesFromLevel(self, clsSet):
        """ Create the SetOfClasses2D from a given iteration. """
        myFileParticles = self._getExtraPath('last_images.xmd')
        myFileClasses = self._getExtraPath('last_classes.xmd')
        self._loadClassesInfo(myFileClasses)
        #blocks = md.getBlocksInMetaDataFile(myFileClasses)
        #for __, block in enumerate(blocks):
        #    if block.startswith('class0'):
        #        xmpMd = block + "@" + myFileClasses
        #        iterator = md.SetMdIterator(xmpMd, sortByLabel=md.MDL_ITEM_ID,
        #                                    updateItemCallback=self._updateParticle,
        #                                    skipDisabled=True)
        #        # itemDataIterator is not necessary because, the class SetMdIterator
        #        # contain all the information about the metadata
        #        clsSet.classifyItems(updateItemCallback=iterator.updateItem,
        #                             updateClassCallback=self._updateClass)
        xmpMd = myFileParticles
        iterator = md.SetMdIterator(xmpMd, sortByLabel=md.MDL_ITEM_ID,
                                    updateItemCallback=self._updateParticle,
                                    skipDisabled=True)

        # itemDataIterator is not necessary because, the class SetMdIterator
        # contain all the information about the metadata
        clsSet.classifyItems(updateItemCallback=iterator.updateItem,
                             updateClassCallback=self._updateClass)


    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        if self.useReferenceImages:
            if self.referenceImages.hasValue():
                refImage = self.referenceImages.get()
                [x1, y1, z1] = refImage.getDim()
                [x2, y2, z2] = self.inputParticles.get().getDim()
                if x1 != x2 or y1 != y2 or z1 != z2:
                    errors.append('The input images and the reference images '
                                  'have different sizes')
            else:
                errors.append("Please, enter the reference images")
        return errors


    def _summary(self):
        summary = []
        if not hasattr(self, 'outputClasses'):
            summary.append("Output alignment not ready yet.")
        else:
            summary.append("Input Particles: %s"
                           % self.inputParticles.get().getSize())
            if self.useReferenceImages:
                summary.append("Aligned with reference images: %s"
                               % self.referenceImages.get().getSize())
            else:
                summary.append("Aligned with no reference images.")
        return summary

    def _citations(self):
        return ['Sorzano2010a']

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputClasses'):
            methods.append("Output alignment not ready yet.")
        else:
            if self.useReferenceImages:
                methods.append(
                    "We aligned images %s with respect to the reference image set "
                    "%s using CL2D [Sorzano2010a]"
                    % (self.getObjectTag('inputParticles'),
                       self.getObjectTag('referenceImages')))
            else:
                methods.append(
                    "We aligned images %s with no references using CL2D "
                    "[Sorzano2010a]" % self.getObjectTag('inputParticles'))
            methods.append(" and produced %s images."
                           % self.getObjectTag('outputClasses'))
        return methods

