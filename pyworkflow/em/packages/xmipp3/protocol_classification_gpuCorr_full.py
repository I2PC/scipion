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
from pyworkflow import VERSION_1_2
from pyworkflow.em import SetOfParticles, ALIGN_2D, ALIGN_NONE
from pyworkflow.em.protocol import ProtAlign2D
import pyworkflow.em.metadata as md
import pyworkflow.protocol.params as params
from pyworkflow.em.metadata.utils import iterRows, getSize
from xmipp import Image, MD_APPEND, DT_DOUBLE
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles, \
    xmippToLocation, rowToAlignment, rowToParticle
from shutil import copy, copytree
from os.path import exists, getmtime
from datetime import datetime
from pyworkflow.utils import prettyTime, cleanPath
from pyworkflow.object import Set
from pyworkflow.protocol.constants import STATUS_NEW
from random import randint
import pyworkflow.protocol.constants as const
from os import system, popen, mkdir, listdir
from os.path import join


HASH_SIZE = 100

class HashTableDict:
    def __init__(self, Ndict=HASH_SIZE):
        self.Ndict = Ndict
        self.dict = [{}]*Ndict

    def isItemPresent(self, idx):
        return idx in self.dict[idx % self.Ndict]

    def pushItem(self, idx):
        idxDict = idx % self.Ndict
        if not idx in self.dict[idxDict]:
            self.dict[idxDict][idx]=1


class XmippProtStrGpuCrrCL2D(ProtAlign2D):
    """ 2D alignment in full streaming using Xmipp GPU Correlation.
    The set of classes will be growing whilst new particle images are
    received."""
    _label = 'gl2d streaming'
    _lastUpdateVersion = VERSION_1_2
    _stepsCheckSecs = 10

    # --------------------------- DEFINE param functions -----------------------
    def _defineAlignParams(self, form):
        form.addParam('maxShift', params.IntParam, default=10,
                      label='Maximum shift (%):',
                      help='Maximum shift allowed during the alignment as '
                           'percentage of the input set size',
                      expertLevel=const.LEVEL_ADVANCED)
        form.addParam('keepBest', params.IntParam, default=2,
                      label='Number of best images:',
                      help='Number of classes to assign every input image '
                           'during the alignment',
                      expertLevel=const.LEVEL_ADVANCED)
        form.addParam('numberOfSplitIterations', params.IntParam, default=2,
                      label='Number of iterations in split stage:',
                      help='Maximum number of iterations in split stage',
                      expertLevel=const.LEVEL_ADVANCED)
        form.addParam('numberOfClassifyIterations', params.IntParam, default=2,
                      label='Number of iterations in classify stage:',
                      help='Maximum number of iterations when the classification'
                           ' of the whole image set is carried out',
                      expertLevel=const.LEVEL_ADVANCED)
        form.addParam('imageSize', params.IntParam, default=64,
                      label='Image size',
                      help='The image size can be downsampled to accelerate '
                           'the classification',
                      expertLevel=const.LEVEL_ADVANCED)
        form.addParam('threshold', params.IntParam, default=500,
                      label='Threshold to split',
                      help='The threshold in the number of images assigned '
                           'to one class to make a spliting of that class',
                      expertLevel=const.LEVEL_ADVANCED)
        form.addParam('blockSize', params.IntParam, default=5000,
                      label='Block size',
                      help='The inputs will be processed in a block-by-block '
                           'basis of this size',
                      expertLevel=const.LEVEL_ADVANCED)
        form.addParam('numClasses', params.IntParam, default=1000,
                      label='Maximum number of classes',
                      help='Maximum number of classes to be generated',
                      expertLevel=const.LEVEL_ADVANCED)
        form.addParallelSection(threads=0, mpi=0)


    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        """" Mainly prepare the command line for calling cuda corrrelation
        program"""
        self._readingCheckPoint()

        self.listInFn = []
        self.listOutFn = []
        self.listOutSplitFn = []
        self.doneListFn = []
        self.lastDate = 0
        self.imgsExp = self._getExtraPath('imagesExp.xmd')
        self.listNumImgs = []
        self.listNameImgs = []
        self.listRefImgs = []
        self.particlesToProcess = []
        self.randRef = None
        self.htAlreadyProcessed = HashTableDict()
        xOrig = self.inputParticles.get().getXDim()
        self.maximumShift = int(self.maxShift.get() * xOrig / 100)

        self._loadInputList()

        self.changeSampling = False
        xO = self.listOfParticles[0].getXDim()
        newSize = self.imageSize.get()
        self.newSamplingRate = self.listOfParticles[0].getSamplingRate()
        if xO != newSize:
            self.changeSampling=True
            factor = float(xO) / float(newSize)
            self.newSamplingRate = self.newSamplingRate * factor

        deps = []
        numBlk, rem = divmod(float(len(self.listOfParticles)),
                             float(self.blockSize.get()))
        numBlk = int(numBlk)
        rem = int(rem)
        if rem > 0:
            numBlk += 1
        for i in range(numBlk):
            if i==0:
                deps += self._insertStepsForParticles(True, False)
            else:
                deps += self._insertStepsForParticles(False, True)

        self._insertFunctionStep('createOutputStep', prerequisities=deps,
                                 wait=True)

    def _insertStepsForParticles(self, flagSplit, reclassification):

        inputImgs = self._getInputFn()
        deps = []

        stepIdClassify = self._insertFunctionStep\
            ('classifyStep', inputImgs, flagSplit, reclassification,
             prerequisites=[])
        deps.append(stepIdClassify)

        return deps


    # --------------------------- STEPS functions --------------------------

    def _stepsCheck(self):
        self._checkNewInput()
        self._checkNewOutput()

    def createOutputStep(self):
        self._createFinalMetadata()
        self._createFinalClasses()

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
        if len(self.listOfParticles)==0:
            return None

        deps=[]
        numBlk, rem = divmod(float(len(self.listOfParticles)),
                             float(self.blockSize.get()))
        numBlk = int(numBlk)
        rem = int(rem)
        if rem>0:
            numBlk+=1
        for i in range(numBlk):
            deps += self._insertStepsForParticles(False, True)

        if outputStep is not None:
            outputStep.addPrerequisites(*deps)
        self.updateSteps()



    def _checkNewOutput(self):
        """ Check for already done files and update the output set. """

        newDone = self._readDoneList()

        # We have finished when there is not more inputs (stream closed)
        # and the number of processed particles is equal to the number of inputs
        self.finished = (self.isStreamClosed == Set.STREAM_CLOSED
                         and len(newDone) == 0)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if newDone:
            self._updateOutputSet(streamMode)
            self.doneListFn += newDone

        elif not self.finished:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)


    def classifyStep(self, expImgMd, flag_split, reclassification):

        auxList=[]
        for im in range(len(self.listOfParticles)):
            idx = self.listOfParticles[im].getObjId()
            if not self.htAlreadyProcessed.isItemPresent(idx):
                auxList.append(self.listOfParticles[im])
                self.htAlreadyProcessed.pushItem(idx)

        self.particlesToProcess+=auxList
        if len(self.particlesToProcess)==0:
            return

        if len(self.particlesToProcess)<self.blockSize.get():
            particlesToProcessAux = self.particlesToProcess
            lastIm = len(self.particlesToProcess)-1
        else:
            particlesToProcessAux = self.particlesToProcess[:self.blockSize.get()]
            lastIm = self.blockSize.get()-1

        self.generateInput(expImgMd, flag_split, reclassification,
                           particlesToProcessAux)

        if flag_split:
            refImgMd = self._getExtraPath('split_last_classes.xmd')
        else:
            refImgMd = self._getExtraPath('last_classes.xmd')

        i=0
        while i <self.numberOfClassifyIterations:
            outImgs, refImgMd = self.iterationStep(refImgMd, expImgMd, i, False)

            if flag_split==False and i+1 < self.numberOfClassifyIterations:
                classesFnPrev = self._getExtraPath('last_classes.xmd')
                self.averageClasses(refImgMd, classesFnPrev, refImgMd, True)
            i += 1

        self.generateOutputClasses(refImgMd, flag_split)

        self.checkSplit()


        self.lastDate = self.particlesToProcess[lastIm].getObjCreation()
        self._saveCreationTimeFile(self.lastDate)

        for p in range(len(particlesToProcessAux)):
            self.particlesToProcess.pop(0)

        self._savingCheckPoint()




    # --------------------------- UTILS functions ------------------------------

    def splitStep(self, expImgMd):
        i=0
        classesOut=None
        while i <self.numberOfSplitIterations:
            outImgs,classesOut = self.iterationStep(classesOut,expImgMd,i,True)
            i+=1
            length = getSize(classesOut)
            if length == 1:
                i = 0

        self.generateMdForClassification(classesOut)



    def generateInput(self, inputImgs, flag_split, reclassification,
                      particlesToProcess):

        inputStack = inputImgs.replace('xmd','stk')
        newSize = self.imageSize.get()
        if self.changeSampling:
            writeSetOfParticles(particlesToProcess, inputImgs,
                                alignType=ALIGN_NONE)
            args = "-i %s -o %s --save_metadata_stack --fourier %d " % (
                inputImgs, inputStack, newSize)
            self.runJob("xmipp_image_resize", args, numberOfMpi=1)
        else:
            writeSetOfParticles(particlesToProcess, inputImgs,
                                alignType=ALIGN_NONE)

        if reclassification and len(self.listNumImgs)!=0:
            self.randRef = randint(1, len(self.listNumImgs))
            self._unionReclassification(self.randRef, inputImgs)

        if flag_split:
            self.splitStep(inputImgs)



    def generateMdForClassification(self, classesOut):

        listNameImgs = self.listNameImgs
        listNumImgs = self.listNumImgs

        count = 1
        # Construct the new classes with the renumerated old classes
        mdNewClasses = md.MetaData()
        for i in range(len(listNumImgs)):
            if listNumImgs[i] != -1:
                name = listNameImgs[i]
                fn = name[name.find('@') + 1:-4] + '.xmd'
                numRef = int(name[0:6])

                mdClass = md.MetaData("classes@" + fn)
                for row in iterRows(mdClass):
                    if mdClass.getValue(md.MDL_REF, row.getObjId()) == numRef:
                        row.setValue(md.MDL_REF, count)
                        row.addToMd(mdNewClasses)
                count += 1

        # Add the two new classes to the list of renumerated classes
        mdClass = md.MetaData("classes@" + classesOut)
        rows = iterRows(mdClass)
        for row in rows:
            row.setValue(md.MDL_REF, count)
            row.addToMd(mdNewClasses)
            count += 1
        mdNewClasses.write('classes@'
                           + self._getExtraPath('split_last_classes.xmd'),
                           MD_APPEND)

        # Generate the intermediate images and the blocks of the intermediate
        # classes for the unchanged classes
        count=1
        for i in range(len(listNumImgs)):
           if listNumImgs[i] != -1:
               # Read the list of images in this class
               mdImgsInClass =  md.MetaData(self._getExtraPath(
                   'dataClass%06d.xmd' % (i+1)))
               mdImgsInClass.fillConstant(md.MDL_REF, count)
               mdImgsInClass.write(self._getExtraPath('dataClass%06d.xmd' %
                                                      count))
               count += 1

        # Add the two new classes
        for newRef in range(0, 2):
           mdImgsInClass = md.MetaData(
               'class%06d_images@%s' % (newRef + 1, classesOut))
           mdImgsInClass.fillConstant(md.MDL_REF, count)
           mdImgsInClass.write(self._getExtraPath('dataClass%06d.xmd' % count))
           count+=1


    def averageClasses(self, finalMetadata, lastMetadata, newMetadata,
                       flag_iter):

        metadataItemLast = md.MetaData(lastMetadata)
        metadataItemNew = md.MetaData(newMetadata)
        if flag_iter:
            metadataFinal = md.MetaData(finalMetadata)
            finalName = metadataFinal.getValue(md.MDL_IMAGE, 1)
            finalName = finalName[7:-3]+'stk'
        else:
            finalName = finalMetadata

        total=[]
        newRef = 1
        i = 1
        for itemLast, itemNew in zip(metadataItemLast, metadataItemNew):
            listToMultiply = []
            numImgsLastClasses = metadataItemLast.getValue(
                md.MDL_CLASS_COUNT, itemLast)
            nameRefLastClasses = metadataItemLast.getValue(md.MDL_IMAGE,
                                                         itemLast)
            nameRefNewClasses = metadataItemNew.getValue(md.MDL_IMAGE,
                                                         itemNew)
            numImgsNewClasses = metadataItemNew.getValue(md.MDL_CLASS_COUNT,
                                                         itemNew)
            if flag_iter==False and i==self.randRef:
                total.append(numImgsNewClasses)
            else:
                total.append(numImgsLastClasses + numImgsNewClasses)
            i += 1

            if (numImgsLastClasses + numImgsNewClasses)==0:
                continue

            if numImgsNewClasses==0:
                im1 = Image(nameRefLastClasses)
                im1.write('%06d@' % newRef + finalName)
                newRef += 1
                continue
            if numImgsLastClasses==0:
                im2 = Image(nameRefNewClasses)
                im2.write('%06d@' % newRef + finalName)
                newRef += 1
                continue

            listToMultiply.append(float(numImgsLastClasses)/
                                  float(numImgsLastClasses + numImgsNewClasses))
            listToMultiply.append(float(numImgsNewClasses)/
                                  float(numImgsLastClasses + numImgsNewClasses))

            im1 = Image(nameRefLastClasses)
            im2 = Image(nameRefNewClasses)
            im1.align(im2)
            im1.inplaceMultiply(listToMultiply[0])
            im2.inplaceMultiply(listToMultiply[1])
            im1.convert2DataType(DT_DOUBLE)
            im2.convert2DataType(DT_DOUBLE)
            im1.inplaceAdd(im2)  # Aligned
            im1.write('%06d@' % newRef + finalName)
            newRef+=1

        return total


    def generateOutputClasses(self, classesOut, firstTime):


        if firstTime:
            self._saveFileDataClasses(classesOut, self._getExtraPath(
               'last_classes.xmd'))
            # Add the two new classes
            for i in range(1, 3):
               mdImgsInClass = md.MetaData(
                   'class%06d_images@%s' % (i,classesOut))
               mdImgsInClass.write(self._getExtraPath('dataClass%06d.xmd' % i))
            return

        finalMetadata = self._getExtraPath('aux_classes.stk')
        lastMetadata = self._getExtraPath('last_classes.xmd')
        newMetadata = classesOut

        total = self.averageClasses(finalMetadata, lastMetadata, newMetadata,
                                    False)

        copy(self._getExtraPath('aux_classes.stk'),
             self._getExtraPath('last_classes.stk'))

        mdAll = md.MetaData()
        newRef=1
        for i in total:
            if i != 0:
                row = md.Row()
                row.setValue(md.MDL_REF, newRef)
                row.setValue(md.MDL_IMAGE, '%06d@'%newRef + finalMetadata)
                row.setValue(md.MDL_CLASS_COUNT, i)
                row.addToMd(mdAll)
                newRef+=1
        mdAll.write('classes@'+finalMetadata[:-3] + 'xmd', MD_APPEND)

        copy(self._getExtraPath('aux_classes.xmd'),
            self._getExtraPath('last_classes.xmd'))

        newRef=1
        for i, val in enumerate(total):
            if val != 0:
                self._unionDataClass(classesOut, i+1, newRef)
                newRef+=1


    def iterationStep (self, refSet, imgsExp, iter, flag_split):

        if flag_split:
            outImgs, classesOut = self._getOutputSplitFn()
        else:
            outImgs, classesOut = self._getOutputClassFn()

        outDirName = imgsExp[:-4]
        if iter==0 and flag_split==True:

            # First step: divide the metadata input file to generate
            # a couple of references
            self._params = {'imgsExp': imgsExp,
                            'outDir': outDirName}
            args = ('-i %(imgsExp)s -n 2 --oroot %(outDir)s')
            self.runJob("xmipp_metadata_split", args % self._params,
                        numberOfMpi=1)

            # Second step: calculate the means of the previous metadata
            expSet1 = outDirName + '000001.xmd'
            avg1 = outDirName + '_000001'
            expSet2 = outDirName + '000002.xmd'
            avg2 = outDirName + '_000002'
            self._params = {'imgsSet': expSet1,
                            'outputAvg': avg1}
            args = ('-i %(imgsSet)s --save_image_stats %(outputAvg)s -v 0')
            self.runJob("xmipp_image_statistics", args % self._params,
                        numberOfMpi=1)

            self._params = {'imgsSet': expSet2,
                            'outputAvg': avg2}
            args = ('-i %(imgsSet)s --save_image_stats %(outputAvg)s -v 0')
            self.runJob("xmipp_image_statistics", args % self._params,
                        numberOfMpi=1)

            # Third step: generate a single metadata with the two previous avgs
            refSet = self._getExtraPath('refSet.xmd')
            self._params = {'avg1': avg1 + 'average.xmp',
                            'avg2': avg2 + 'average.xmp',
                            'outputMd': refSet}
            args = ('-i %(avg1)s --set union %(avg2)s -o %(outputMd)s')
            self.runJob("xmipp_metadata_utilities", args % self._params,
                        numberOfMpi=1)

        # Fourth step: calling program xmipp_cuda_correlation
        if flag_split:
            self._params = {'imgsRef': refSet,
                            'imgsExp': imgsExp,
                            'maxshift': self.maximumShift,
                            'Nrefs': getSize(refSet),
                            'outDir': self._getExtraPath(),
                            'rootFn': classesOut.split('/')[-1].replace(
                                '.xmd','')
                            }
            args = '-i %(imgsExp)s --ref0 %(imgsRef)s --nref %(Nrefs)d ' \
                   '--iter 1 --distance correlation --classicalMultiref ' \
                   '--maxShift %(maxshift)d --odir %(outDir)s --oroot %(' \
                   'rootFn)s'

            self.runJob("mpirun -np 4 -bynode xmipp_mpi_classify_CL2D",
                        args % self._params)

            fileTocopy = classesOut.replace('.xmd','_classes.xmd')
            fileTocopy = fileTocopy.replace('extra/', 'extra/level_00/')
            copy(fileTocopy, classesOut)
            copy(self._getExtraPath("images.xmd"), outImgs)
        else:
            self._params = {'imgsRef': refSet,
                            'imgsExp': imgsExp,
                            'outputFile': outImgs,
                            'keepBest': self.keepBest.get(),
                            'maxshift': self.maximumShift,
                            'outputClassesFile': classesOut
                            }
            args = '-i_ref %(imgsRef)s -i_exp %(imgsExp)s -o %(outputFile)s '\
                   '--keep_best %(keepBest)d --maxShift %(maxshift)d ' \
                   '--classify %(outputClassesFile)s --simplifiedMd'
            self.runJob("xmipp_cuda_correlation", args % self._params,
                        numberOfMpi=1)

        if exists(outDirName + '000001.xmd'):
            cleanPath(expSet1)
            cleanPath(expSet2)
            cleanPath(avg1 + 'average.xmp')
            cleanPath(avg2 + 'average.xmp')
            cleanPath(avg1 + 'stddev.xmp')
            cleanPath(avg2 + 'stddev.xmp')

        return outImgs, classesOut


    def checkSplit(self):

        outSet = self._getExtraPath('last_classes.xmd')

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

        i=0
        auxList = sorted(self.listNumImgs, reverse=True)
        while i<len(self.listNumImgs):

            if len(self.listNumImgs) < self.numClasses.get(): #inside the while
                # just in case we lose some class we want to allow one more
                # split
                if auxList[i] < (1.75 * self.threshold.get()):
                    i+=1
                    continue

                maxValue = auxList[i]
                maxPos = self.listNumImgs.index(maxValue)

                self.listNumImgs[maxPos] = -1
                bestRef = self.listRefImgs[maxPos]

                outputMd = self._getExtraPath('dataClass%06d.xmd' % bestRef)
                self.splitStep(outputMd)

                i=0
                outSet = self._getExtraPath('split_last_classes.xmd')

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

                copy(outSet, self._getExtraPath('last_classes.xmd'))
                auxList = sorted(self.listNumImgs, reverse=True)
                i=0

            else:
                break


    def _loadInputList(self):
        """ Load the input set of ctfs and create a list. """
        particlesSet = self._loadInputParticleSet()
        lastDate = self._readCreationTimeFile()
        self.isStreamClosed = particlesSet.getStreamState()
        self.listOfParticles = []

        for p in particlesSet.iterItems(orderBy='creation',
                                        where="creation>'%s'"
                                        % lastDate):
            idx = p.getObjId()
            if not self.htAlreadyProcessed.isItemPresent(idx):
                newPart = p.clone()
                newPart.setObjCreation(p.getObjCreation())
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

    def _getOutputSplitFn(self):
        nameImages = self._getExtraPath('split_general_images')
        imagesFn = self._getUniqueFn(nameImages, self.listOutSplitFn)
        classesFn = imagesFn.replace('images', 'classes')
        return imagesFn, classesFn

    def _getOutputClassFn(self):
        nameImages = self._getExtraPath('class_general_images')
        imagesFn = self._getUniqueFn(nameImages, self.listOutFn)
        classesFn = imagesFn.replace('images', 'classes')
        return imagesFn, classesFn

    def _readDoneList(self):
        return [fn for fn in self.listOutFn if fn not in self.doneListFn]

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == 'createOutputStep':
                return s
        return None


    def _createFinalClasses(self):

        inputs = self.inputParticles.get()
        # Here the defineOutputs function will call the write() method
        outputSet = self._createSetOfClasses2D(inputs)
        self._fillClasses(outputSet)
        self._defineOutputs(**{'outputClasses': outputSet})
        self._defineSourceRelation(self.inputParticles, outputSet)
        self._store(outputSet)
        # Close set databaset to avoid locking it
        outputSet.close()

    def _updateOutputSet(self, state=Set.STREAM_OPEN):

        setFile = self._getPath('averages.sqlite')
        if exists(setFile):
            cleanPath(setFile)
        inputs = self.inputParticles.get()
        # Here the defineOutputs function will call the write() method
        outputAvg = self._createSetOfAverages()
        outputAvg.copyInfo(inputs)
        self._fillAverages(outputAvg)
        self._defineOutputs(**{'outputAverages': outputAvg})
        self._defineSourceRelation(self.inputParticles, outputAvg)
        self._store(outputAvg)
        # Close set databaset to avoid locking it
        outputAvg.close()


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
            rep.setSamplingRate(self.newSamplingRate)

    def _loadClassesInfo(self, filename):
        """ Read some information about the produced 2D classes
        from the metadata file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id
        mdClasses = md.MetaData(filename)
        for classNumber, row in enumerate(md.iterRows(mdClasses)):
            index, fn = xmippToLocation(row.getValue(md.MDL_IMAGE))
            self._classesInfo[classNumber + 1] = (index, fn, row.clone())


    def _fillClasses(self, clsSet):
        """ Create the SetOfClasses2D from a given iteration. """
        myFileParticles = self._getExtraPath('last_images.xmd')
        myFileClasses = self._getExtraPath('last_classes.xmd')
        self._loadClassesInfo(myFileClasses)
        xmpMd = myFileParticles
        iterator = md.SetMdIterator(xmpMd, sortByLabel=md.MDL_ITEM_ID,
                                    updateItemCallback=self._updateParticle,
                                    skipDisabled=True)
        clsSet.classifyItems(updateItemCallback=iterator.updateItem,
                             updateClassCallback=self._updateClass)


    def _fillAverages(self, avgSet):
        """ Create the SetOfAverages from a given metadata """
        myFileClasses = "classes@" + self._getExtraPath('last_classes.xmd')
        repSet = md.MetaData(myFileClasses)
        for rep in iterRows(repSet):
            particle = rowToParticle(rep)
            repId = rep.getValue(md.MDL_REF) #rep.getObjId()
            particle.setObjId(repId)
            avgSet.append(particle)


    def _saveFileDataClasses(self, fn, fnSave):

        line = ""
        numLine = 0
        while not line.startswith('data_class0'):
            numLine += 1
            line = popen('sed -n %ip %s' % (numLine, fn)).read()
        system('head -%i %s >> %s' % (numLine - 1, fn, fnSave))


    def _createFinalMetadata(self):
        fnClasses = self._getExtraPath('last_classes.xmd')
        numClasses = md.getSize(fnClasses)

        for i in range(numClasses):
            num=i+1
            fn = self._getExtraPath('dataClass%06d.xmd' % num)
            if exists(fn):
                line = ""
                numLine = 0
                while not line.startswith('    '):
                    numLine += 1
                    line = popen('sed -n %ip %s' % (numLine, fn)).read()
                aux = self._getExtraPath('aux.xmd')
                system('head -%i %s >> %s' % (numLine - 1, fn, aux))

                dataHeader = self._getExtraPath('dataHeader.xmd')
                fp = open(aux, 'r')
                fout = open(dataHeader, 'w')
                for i, line in enumerate(fp):
                    if i<2:
                        continue
                    if not line.startswith('data_noname'):
                        fout.write(line)
                    else:
                        line = line.replace('noname', 'class%06d_images'%num)
                        fout.write(line)
                fp.close()
                fout.close()
                cleanPath(aux)
                fnOut = self._getExtraPath('dataImages%06d.xmd' % num)
                system('tail -n +%i %s >> %s' % (i+2, fn, aux))
                system('cat %s %s >> %s'%(dataHeader, aux, fnOut))
                if num==1:
                    system('cat %s >> %s' % (fn,
                                        self._getExtraPath('last_images.xmd')))
                else:
                    system('cat %s >> %s' %(aux,
                                        self._getExtraPath('last_images.xmd')))
                cleanPath(aux)
                cleanPath(dataHeader)
                system('cat %s >> %s' % (fnOut,
                                        self._getExtraPath('last_classes.xmd')))


    def _unionDataClass(self, fn, refOld, refNew):

        fnNew = self._getExtraPath('dataClass%06d.xmd' % refNew)
        fnOld = self._getExtraPath('dataClass%06d.xmd' % refOld)
        if refOld != refNew and exists(fnOld):
            mdOld = md.MetaData(fnOld)
            mdOld.fillConstant(md.MDL_REF, refNew)
            mdOld.write(fnOld)

        if exists(fnOld):
            fnAux = self._getExtraPath('aux.xmd')
            mdImgsInClass = md.MetaData('class%06d_images@%s' % (refOld, fn))
            if not mdImgsInClass.isEmpty():
                if refOld != refNew:
                    mdImgsInClass.fillConstant(md.MDL_REF, refNew)
                mdImgsInClass.write(fnAux)
                line = ""
                numLine = 0
                while not line.startswith('   '):
                    numLine += 1
                    line = popen('sed -n %ip %s' % (numLine, fnAux)).read()
                fnSave = self._getExtraPath('newDataClass%06d.xmd'%refNew)
                system('tail -n +%i %s >> %s' % (numLine, fnAux, fnSave))
                cleanPath(fnAux)
                if refOld==refNew:
                    system('cat %s >> %s' % (fnSave, fnNew))
                else:
                    if exists(fnNew):
                        cleanPath(fnNew)
                    system('cat %s %s >> %s' %(fnOld, fnSave, fnNew))
                cleanPath(fnSave)
            else:
                if refOld!=refNew:
                    if exists(fnNew):
                        cleanPath(fnNew)
                    copy(fnOld, fnNew)
        else:
            mdImgsInClass = md.MetaData('class%06d_images@%s' % (refOld, fn))
            if refOld != refNew:
                mdImgsInClass.fillConstant(md.MDL_REF, refNew)
            mdImgsInClass.write(fnNew)


    def _unionReclassification(self, ref, inputImgs):

        fn1 = self._getExtraPath('dataClass%06d.xmd' % ref)
        line = ""
        numLine = 0
        while not line.startswith('   '):
            numLine += 1
            line = popen('sed -n %ip %s' % (numLine, fn1)).read()
        fnSave = self._getExtraPath('aux.xmd')
        system('tail -n +%i %s >> %s' % (numLine, fn1, fnSave))
        system('cat %s >> %s' % (fnSave, inputImgs))
        cleanPath(fnSave)
        cleanPath(fn1)


    def _saveCreationTimeFile(self, cTime):
        fn = open(self._getExtraPath('creation.txt'),'w')
        fn.write(cTime)
        fn.close()

    def _readCreationTimeFile(self):
        if exists(self._getExtraPath('creation.txt')):
            fn = open(self._getExtraPath('creation.txt'), 'r')
            cTime = fn.readline()
            fn.close()
        else:
            cTime = 0
        return cTime

    def _savingCheckPoint(self):
        if not exists(join(self._getExtraPath(), 'checkpoint')):
            mkdir(join(self._getExtraPath(), 'checkpoint'))

        listFolder = listdir(self._getExtraPath())
        for fn in listFolder:
            if fn.startswith('dataClass'):
                copy(join(self._getExtraPath(),fn),
                     self._getExtraPath(join('checkpoint',fn)))
        copy(self._getExtraPath('last_classes.xmd'),
             self._getExtraPath(join('checkpoint', 'last_classes.xmd')))
        if exists(self._getExtraPath('last_classes.stk')):
            copy(self._getExtraPath('last_classes.stk'),
                 self._getExtraPath(join('checkpoint', 'last_classes.stk')))

    def _readingCheckPoint(self):
        if exists(join(self._getExtraPath(), 'checkpoint')):
            listFolder = listdir(join(self._getExtraPath(), 'checkpoint'))
            for fn in listFolder:
                copy(self._getExtraPath(join('checkpoint',fn)),
                     self._getExtraPath(fn))




    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        newSize = self.imageSize.get()
        x, y, _ = self.inputParticles.get().getDim()
        if newSize>x or newSize>y:
            errors.append('The image size must be smaller than the size of '
                          'the input images')
        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputClasses'):
            summary.append("Output alignment not ready yet.")
        else:
            summary.append("Input Particles: %s"
                           % self.inputParticles.get().getSize())
        summary.append("Aligned in streaming.")
        return summary

    def _citations(self):
        return ['Sorzano2010a']

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputClasses'):
            methods.append("Output alignment not ready yet.")
        else:
            methods.append(
                "We aligned images %s in streaming using CL2D "
                "[Sorzano2010a] and GPU correlation methods " %
                self.getObjectTag('inputParticles'))
            methods.append(" and produced %s images."
                           % self.getObjectTag('outputClasses'))
        return methods

