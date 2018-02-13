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
from pyworkflow.em import SetOfParticles, SetOfClasses2D, ALIGN_2D, ALIGN_NONE
from pyworkflow.em.protocol import ProtAlign2D
import pyworkflow.em.metadata as md
import pyworkflow.protocol.params as params
from pyworkflow.em.metadata.utils import iterRows, getSize
from xmipp import Image, MD_APPEND, DT_DOUBLE
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles, \
    xmippToLocation, rowToAlignment, rowFromMd, rowToClass, rowToParticle
from shutil import copy
from os.path import exists, getmtime
from datetime import datetime
from pyworkflow.utils import prettyTime, cleanPath
from pyworkflow.object import Set
from pyworkflow.protocol.constants import STATUS_NEW
from random import randint
import time
from pyworkflow.monitor import Timer


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
    """ Aligns a set of particles in streaming using the GPU Correlation
    algorithm. """
    _label = 'align with GPU Correlation in full streaming'
    _lastUpdateVersion = VERSION_1_2
    _stepsCheckSecs = 10

    # --------------------------- DEFINE param functions -----------------------
    def _defineAlignParams(self, form):
        form.addParam('maximumShift', params.IntParam, default=10,
                      label='Maximum shift (px):')
        form.addParam('keepBest', params.IntParam, default=2,
                      label='Number of best images:',
                      help='Number of the best images to keep for every class')
        form.addParam('numberOfSplitIterations', params.IntParam, default=1,
                      label='Number of iterations in split stage:',
                      help='Maximum number of iterations in split stage')
        form.addParam('numberOfClassifyIterations', params.IntParam, default=1,
                      label='Number of iterations in classify stage:',
                      help='Maximum number of iterations when the classification'
                           ' of the whole image set is carried out')
        form.addParam('imageSize', params.IntParam, default=64,
                      label='Image size',
                      help='The image size can be downsampled to accelerate '
                           'the classification')
        form.addParallelSection(threads=0, mpi=0)


    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        """" Mainly prepare the command line for calling cuda corrrelation program"""

        self.listInFn = []
        self.listOutFn = []
        self.listOutSplitFn = []
        self.doneListFn = []
        self.lastDate = 0
        self.imgsExp = self._getExtraPath('imagesExp.xmd')
        self.listNumImgs = []
        self.listNameImgs = []
        self.listRefImgs = []
        self.htAlreadyProcessed = HashTableDict()
        self.last_time = time.time()
        self.deltaTime = 60.0

        self._loadInputList()
        deps = []
        deps += self._insertStepsForParticles(True, False)
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
        self._checkNewOutput()

    def _checkNewInput(self):
        """ Check if there are new particles to be processed and add
        the necessary steps."""
        initTime = time.time()
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
        deps += self._insertStepsForParticles(False, True)
        if outputStep is not None:
            outputStep.addPrerequisites(*deps)
        self.updateSteps()

        final_time = time.time()
        exec_time = final_time-initTime
        print("_checkNewInput exec_time", exec_time)


    def _checkNewOutput(self):
        """ Check for already done files and update the output set. """
        inputSize = self.inputParticles.get().getSize()
        if inputSize>10000 and inputSize<50000:
            self.deltaTime = 180.0
        elif inputSize>=50000 and inputSize<100000:
            self.deltaTime = 300.0
        elif inputSize>100000:
            self.deltaTime = 600.0

        initial_time = time.time()
        if initial_time < self.last_time + self.deltaTime:
            return
        else:
            self.last_time = initial_time

        print("En checkNewOutput")

        newDone = self._readDoneList()

        # We have finished when there is not more inputs (stream closed)
        # and the number of processed particles is equal to the number of inputs
        self.finished = (self.isStreamClosed == Set.STREAM_CLOSED
                         and len(newDone) == 0)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if newDone:
            print("Hay newDone")
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

        final_time = time.time()
        exec_time = final_time - initial_time
        print("_checkNewOutput exec_time", exec_time)


    def classifyStep(self, expImgMd, flag_split, reclassification):

        print("len(self.listOfParticles)", len(self.listOfParticles))

        self.generateInput(expImgMd, flag_split, reclassification)

        for p in self.listOfParticles:
            partId = p.getObjId()
            self.htAlreadyProcessed.pushItem(partId)

        if flag_split:
            refImgMd = self._getExtraPath('split_last_classes.xmd')
        else:
            refImgMd = self._getExtraPath('last_classes.xmd')

        i=0
        while i <self.numberOfClassifyIterations:
            print("ITER",i)
            outImgs, refImgMd = self.iterationStep(refImgMd, expImgMd, i, False)

            if flag_split==False and i+1 < self.numberOfClassifyIterations:
                classesFnPrev = self._getExtraPath('last_classes.xmd')
                self.averageClasses(refImgMd, classesFnPrev, refImgMd, True)
            i += 1

        self.generateOutputClasses(refImgMd, outImgs, flag_split)
        self.checkSplit(refImgMd)



    # --------------------------- UTILS functions ------------------------------

    def splitStep(self, expImgMd):

        print("len(self.listOfParticles)",len(self.listOfParticles))

        i=0
        clasesOut=None
        while i <self.numberOfSplitIterations:
            outImgs,clasesOut = self.iterationStep(clasesOut,expImgMd,i,True)
            i+=1
            length = getSize(clasesOut)
            if length == 1:
                i = 0
        self.generateMdForClassification(self.listNameImgs, self.listNumImgs,
                                         clasesOut)


    def generateInput(self, inputImgs, flag_split, reclassification):

        inputStack = inputImgs.replace('xmd','stk')
        xO = self.listOfParticles[0].getXDim()
        newSize = self.imageSize.get()
        self.newSamplingRate = self.listOfParticles[0].getSamplingRate()
        if xO != newSize:
            factor = float(xO) / float(newSize)
            self.newSamplingRate = self.newSamplingRate * factor
            writeSetOfParticles(self.listOfParticles, inputImgs,
                                alignType=ALIGN_NONE)
            args = "-i %s -o %s --save_metadata_stack --fourier %d " % (
                inputImgs, inputStack, newSize)
            self.runJob("xmipp_image_resize", args, numberOfMpi=1)
        else:
            writeSetOfParticles(self.listOfParticles, inputImgs,
                                alignType=ALIGN_NONE)

        if reclassification:
            randRef = randint(1, len(self.listNumImgs))
            outMd = self._getExtraPath('last_classes.xmd')
            block = 'class%06d' % (randRef)
            self._params = {'newMd': block + "_images@" + outMd,
                            'outMd': inputImgs,
                            }
            args = ('-i %(newMd)s -o %(outMd)s --set union_all %(outMd)s')
            self.runJob("xmipp_metadata_utilities",
                        args % self._params, numberOfMpi=1)

            self._params = {'newMd': self._getExtraPath('last_images.xmd'),
                            'label': randRef,
                            }
            args=('-i %(newMd)s -o %(newMd)s --query select "ref != %(label)d"')
            self.runJob("xmipp_metadata_utilities",
                        args % self._params, numberOfMpi=1)

            mdAll = md.MetaData()
            fn = self._getExtraPath('last_classes.xmd')
            mdClass = md.MetaData('classes' + "@" + fn)
            for row in iterRows(mdClass):
                if mdClass.getValue(md.MDL_REF, row.getObjId())==randRef:
                    row.setValue(md.MDL_CLASS_COUNT, 0L)
                row.addToMd(mdAll)
            mdAll.write('classes@' + fn, MD_APPEND)

            count = 1
            blocks = md.getBlocksInMetaDataFile(fn)
            for block in blocks:
                mdAll2 = md.MetaData()
                if block.startswith('class00'):
                    if not block.startswith('class%06d' % (randRef)):
                        mdClass = md.MetaData(block + "@" + fn)
                        rows = iterRows(mdClass)
                        for row in rows:
                            row.addToMd(mdAll2)
                    mdAll2.write('class%06d' % (count) + '_images@' + fn,
                                 MD_APPEND)
                    count+=1

        if flag_split:
            self.splitStep(inputImgs)


    def generateMdForClassification(self, listNameImgs, listNumImgs, clasesOut):

        # Renumerate unchanged classes
        listNewNumImages=[-1]*len(listNumImgs)
        count = 1
        for i in range(len(listNumImgs)):
            if listNumImgs[i] is not -1:
                listNewNumImages[i] = count
                count+=1

        # Construct the new classes with the renumerated old classes
        mdNewClasses = md.MetaData()
        for i in range(len(listNumImgs)):
            if listNumImgs[i] is not -1:
                name = listNameImgs[i]
                #numRef = listRefImgs[i]
                fn = name[name.find('@') + 1:-4] + '.xmd'
                numRef = int(name[0:6])
                print(name, numRef)

                mdClass = md.MetaData("classes@" + fn)
                for row in iterRows(mdClass):
                    if mdClass.getValue(md.MDL_REF, row.getObjId()) == numRef:
                        row.setValue(md.MDL_REF, listNewNumImages[i])
                        row.addToMd(mdNewClasses)

        # Add the two new classes to the list of renumerated classes
        mdClass = md.MetaData("classes@" + clasesOut)
        rows = iterRows(mdClass)
        for row in rows:
            row.setValue(md.MDL_REF, count)
            row.addToMd(mdNewClasses)
            count = count + 1
        mdNewClasses.write('classes@'
                           + self._getExtraPath('split_last_classes.xmd'),
                           MD_APPEND)

        # Generate the intermediate images and the blocks of the intermediate
        # classes for the unchanged classes
        mdAll = md.MetaData()
        for i in range(len(listNumImgs)):
            if listNumImgs[i] is not -1:
                name = listNameImgs[i]
                #numRef = listRefImgs[i]
                fn = name[name.find('@')+1:-4]+'.xmd'
                numRef = int(name[0:6])
                print(name, numRef)

                # Read the list of images in this class
                mdImgsInClass = md.MetaData('class%06d_images@%s' % (numRef,fn))
                mdImgsInClass.fillConstant(md.MDL_REF,listNewNumImages[i])
                mdImgsInClass.write('class%06d' % (listNewNumImages[i]) +
                                    '_images@' + self._getExtraPath(
                                    'split_last_classes.xmd'), MD_APPEND)
                mdAll.unionAll(mdImgsInClass)

        # Add the two new classes
        if len(listNumImgs)==0:
            count=1
        else:
            count=len(listNumImgs)
        for newRef in range(0,2):
            mdImgsInClass = md.MetaData('class%06d_images@%s' % (newRef+1,
                                                                 clasesOut))
            mdImgsInClass.fillConstant(md.MDL_REF,count)
            mdImgsInClass.write('class%06d' % (count) + '_images@' +
                                self._getExtraPath('split_last_classes.xmd'),
                                MD_APPEND)
            mdAll.unionAll(mdImgsInClass)
            count = count+1

        # Write the list of images with their new reference
        mdAll.write(self._getExtraPath('split_last_images.xmd'))


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
        i = 0
        for item in metadataItemLast:
            listToMultiply = []
            numImgsLastClasses = metadataItemLast.getValue(
                md.MDL_CLASS_COUNT, item)
            nameRefLastClasses = metadataItemLast.getValue(md.MDL_IMAGE, item)
            nameRefNewClasses = metadataItemNew.getValue(md.MDL_IMAGE, item)
            numImgsNewClasses = metadataItemNew.getValue(md.MDL_CLASS_COUNT,
                                                        item)
            total.append(numImgsLastClasses + numImgsNewClasses)

            print("numImgsNewClasses", numImgsNewClasses)
            print("numImgsLastClasses", numImgsLastClasses)

            if numImgsNewClasses==0:
                im1 = Image(nameRefLastClasses)
                im1.write('%06d@' % (i + 1) + finalName)
                i += 1
                continue
            if numImgsLastClasses==0:
                im2 = Image(nameRefNewClasses)
                im2.write('%06d@' % (i + 1) + finalName)
                i += 1
                continue

            if total[i]==0: #AJ to avoid dividing by zero
                listToMultiply = [0, 0]
            else:
                listToMultiply.append(float(numImgsLastClasses)/float(total[i]))
                listToMultiply.append(float(numImgsNewClasses)/float(total[i]))

            im1 = Image(nameRefLastClasses)
            im2 = Image(nameRefNewClasses)
            im1.inplaceMultiply(listToMultiply[0])
            im2.inplaceMultiply(listToMultiply[1])
            im1.align(im2)
            im1.convert2DataType(DT_DOUBLE)
            im2.convert2DataType(DT_DOUBLE)
            im1.inplaceAdd(im2)  # Aligned
            im1.write('%06d@' % (i + 1) + finalName)
            i+=1

        return total


    def generateOutputClasses(self, clasesOut, outImgs, firstTime):

        print("En generateOutputClasses")

        if firstTime:
            copy(clasesOut, self._getExtraPath('last_classes.xmd'))
            copy(outImgs, self._getExtraPath('last_images.xmd'))
            return

        finalMetadata = self._getExtraPath('aux_classes.stk')
        lastMetadata = self._getExtraPath('last_classes.xmd')
        newMetadata = clasesOut

        total = self.averageClasses(finalMetadata, lastMetadata, newMetadata,
                                    False)

        copy(self._getExtraPath('aux_classes.stk'),
             self._getExtraPath('last_classes.stk'))

        mdAll = md.MetaData()
        for i in range(len(total)):
            row = md.Row()
            row.setValue(md.MDL_REF, i+1)
            row.setValue(md.MDL_IMAGE, '%06d@'%(i+1) + finalMetadata)
            row.setValue(md.MDL_CLASS_COUNT, total[i])
            row.addToMd(mdAll)
        mdAll.write('classes@'+finalMetadata[:-3] + 'xmd', MD_APPEND)

        for i in range(len(total)):
            self._params = {'lastMd': 'class%06d_images@' % (i + 1) +
                                      lastMetadata,
                            'newMd': 'class%06d_images@' % (i + 1) +
                                     newMetadata,
                            'outMd': 'class%06d_images@' % (i + 1) +
                                     finalMetadata[:-3] + 'xmd'}
            args = ('-i %(lastMd)s --set union_all %(newMd)s '
                    '-o %(outMd)s --mode append')
            self.runJob("xmipp_metadata_utilities",
                        args % self._params, numberOfMpi=1)

        copy(self._getExtraPath('aux_classes.xmd'),
             self._getExtraPath('last_classes.xmd'))

        self._params = {'lastMd': self._getExtraPath('last_images.xmd'),
                        'newMd': outImgs,
                        'outMd': self._getExtraPath('aux_images.xmd')}
        args = ('-i %(lastMd)s --set union %(newMd)s -o %(outMd)s')
        self.runJob("xmipp_metadata_utilities",
                    args % self._params, numberOfMpi=1)

        copy(self._getExtraPath('aux_images.xmd'),
             self._getExtraPath('last_images.xmd'))


    def iterationStep (self, refSet, imgsExp, iter, flag_split):

        if flag_split:
            outImgs, clasesOut = self._getOutputSplitFn()
        else:
            outImgs, clasesOut = self._getOutputClassFn()

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
                            'maxshift': self.maximumShift.get(),
                            'Nrefs': getSize(refSet),
                            'outDir': self._getExtraPath(),
                            'rootFn': clasesOut.split('/')[-1].replace(
                                '.xmd','')
                            }
            args = '-i %(imgsExp)s --ref0 %(imgsRef)s --nref %(Nrefs)d ' \
                   '--iter 1 --distance correlation --classicalMultiref ' \
                   '--maxShift %(maxshift)d --odir %(outDir)s --oroot %(' \
                   'rootFn)s'

            self.runJob("mpirun -np 4 -bynode xmipp_mpi_classify_CL2D",
                        args % self._params)

            fileTocopy = clasesOut.replace('.xmd','_classes.xmd')
            fileTocopy = fileTocopy.replace('extra/', 'extra/level_00/')
            copy(fileTocopy, clasesOut)
            copy(self._getExtraPath("images.xmd"), outImgs)
        else:
            self._params = {'imgsRef': refSet,
                            'imgsExp': imgsExp,
                            'outputFile': outImgs,
                            'keepBest': self.keepBest.get(),
                            'maxshift': self.maximumShift.get(),
                            'outputClassesFile': clasesOut
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

        return outImgs, clasesOut


    def checkSplit(self, refImgMd):

        print("En checkSplit")

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

        thGlobal = 500
        i=0
        while i<len(self.listNumImgs):

            print("self.listNameImgs", self.listNameImgs)
            print("self.listNumImgs", self.listNumImgs)
            print("self.listRefImgs", self.listRefImgs)

            if self.listNumImgs[i]<1.75*thGlobal:
                i+=1
                continue

            maxValue = self.listNumImgs[i]
            maxPos = self.listNumImgs.index(maxValue)

            self.listNumImgs[maxPos] = -1
            bestRef = self.listRefImgs[maxPos]

            outputMd = refImgMd.replace('.xmd', '_major.xmd')

            self._params = {'input': 'class%06d_images' % (bestRef)+'@'+outSet,
                            'outputMd': outputMd
                            }
            args = ('-i %(input)s -o %(outputMd)s')
            self.runJob("xmipp_metadata_utilities", args %
                        self._params, numberOfMpi=1)

            #Split with XXXX_major.xmd
            #Generate metadata with the output of the previous split and the
            # data in general_XXX_classes.xmd (removing the major class)
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
            copy(self._getExtraPath('split_last_images.xmd'),
                 self._getExtraPath('last_images.xmd'))
            refImgMd = self._getExtraPath('last_classes.xmd')

        print("self.listNameImgs", self.listNameImgs)
        print("self.listNumImgs", self.listNumImgs)
        print("self.listRefImgs", self.listRefImgs)


    def _loadInputList(self):
        """ Load the input set of ctfs and create a list. """
        with Timer() as t:
            particlesSet = self._loadInputParticleSet()
        print "_loadInputParticleSet: %s s" % t.secs

        self.isStreamClosed = particlesSet.getStreamState()
        self.listOfParticles = []
        with Timer() as t:
            for p in particlesSet.iterItems(orderBy='creation',
                                            where="creation>'%s'"
                                            % self.lastDate):
                idx = p.getObjId()
                if not self.htAlreadyProcessed.isItemPresent(idx):
                    newPart = p.clone()
                    self.listOfParticles.append(newPart)
            if len(self.listOfParticles)>0:
                self.lastDate = p.getObjCreation()
        print "_loop: %s s" % t.secs

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


    def _updateOutputSet(self, state=Set.STREAM_OPEN):

        setFile = self._getPath('classes2D.sqlite')
        if exists(setFile):
            cleanPath(setFile)

        inputs = self.inputParticles.get()
        # Here the defineOutputs function will call the write() method
        outputSet = self._createSetOfClasses2D(inputs)
        self._fillClasses(outputSet)
        self._defineOutputs(**{'outputClasses': outputSet})
        self._defineSourceRelation(self.inputParticles, outputSet)
        self._store(outputSet)

        # Close set databaset to avoid locking it
        outputSet.close()

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

