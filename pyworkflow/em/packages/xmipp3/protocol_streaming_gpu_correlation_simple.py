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

from pyworkflow.em import SetOfParticles, SetOfClasses2D, ALIGN_2D
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
        form.addParam('numberOfClassifyIterations', params.IntParam, default=1,
                      label='Number of iterations in classify stage:',
                      help='Maximum number of iterations when the classification '
                           'of the whole image set is carried out')
        form.addParallelSection(threads=0, mpi=4)


    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        """" Mainly prepare the command line for calling cuda corrrelation program"""

        self.imgsExp = self._getExtraPath('imagesExp.xmd')
        self.imgsRef = self._getExtraPath('imagesRef.xmd')
        self.level=0
        self.blockToProc=300
        self.listToProc=[]
        self.totalInList=0
        self.numberProcessed = 0
        self.htAlreadyProcessed = HashTableDict()

        self.last_time = time.time()

        self._loadInputList()
        stepsToAppend = self.divideInBlocks()

        fDeps=[]
        self._insertFunctionStep('_convertSetStep', self.imgsRef)
        for i in range(stepsToAppend):
            fDeps += self._insertStepsForParticles()

        self._insertFunctionStep('createOutputStep', prerequisities=fDeps, wait=True)


    def _insertStepsForParticles(self):
        deps = []

        stepIdClassify = self._insertFunctionStep('_classifyStep')
        deps.append(stepIdClassify)

        return deps


    # --------------------------- STEPS functions --------------------------

    def _stepsCheck(self):
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        """ Check if there are new particles to be processed and add the necessary
        steps."""
        initial_time = time.time()
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

        if len(self.listOfParticles) > 0:
            stepsToAppend = self.divideInBlocks()
            fDeps=[]
            init2 = time.time()
            for i in range(stepsToAppend):
                fDeps += self._insertStepsForParticles()
            finish2 = time.time()
            init3 = time.time()
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()
            finish3=time.time()
            print("Calculo 1 exec time", finish2-init2)
            print("Calculo 2 exec time", finish3 - init3)

        final_time = time.time()
        exec_time = final_time-initial_time
        print("_checkNewInput exec_time", exec_time)


    def _checkNewOutput(self):
        """ Check for already done files and update the output set. """
        # Load previously done items (from text file)

        initial_time = time.time()
        #particlesListId = self._readParticlesId()
        particlesListId = self.numberProcessed

        #if particlesListId<100000:
        #    interval = 30.0#180.0
        #else:
        #    interval = 300.0
        #if initial_time<self.last_time+interval:
        #    return
        #else:
        #    self.last_time = initial_time
        #    print("last_time",self.last_time)

        doneList = self._readDoneList()

        # Check for newly done items
        if len(doneList)==0:
            newDone = range(0, particlesListId + 1)
        else:
            newDone = range(doneList[0], particlesListId+1)

        # We have finished when there is not more inputs (stream closed)
        # and the number of processed particles is equal to the number of inputs
        self.finished = (self.isStreamClosed == Set.STREAM_CLOSED
                         and len(newDone)==0)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if len(newDone)>1:
            self._writeDoneList(newDone[len(newDone)-1])
        elif not self.finished:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        outSet = self._loadOutputSet(SetOfClasses2D, 'classes2D.sqlite')
        self._updateOutputSetOfClasses('outputClasses', outSet, streamMode)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

        if (exists(self._getPath('classes2D.sqlite'))):
            outSet.close()

        if exists(self._getExtraPath('last_images.xmd')):
            remove(self._getExtraPath('last_images.xmd'))
        if exists(self._getExtraPath('last_classes.xmd')):
            remove(self._getExtraPath('last_classes.xmd'))

        final_time = time.time()
        exec_time = final_time - initial_time
        print("_checkNewOutput exec_time", exec_time)


    def _convertSetStep(self, imgsFileNameClasses):

        if self.useAsRef==REF_CLASSES:
            setOfClasses = self.inputClasses
        else:
            setOfClasses = self.inputAverages

        self.num_classes=setOfClasses.get().__len__()
        if self.useAsRef == REF_CLASSES:
            writeSetOfClasses2D(setOfClasses.get(),
                                imgsFileNameClasses,
                                writeParticles=True)
        else:
            writeSetOfParticles(setOfClasses.get(), imgsFileNameClasses)


    def _classifyStep(self):

        #time.sleep(30)

        self._generateInputMd()

        expImgMd = self._getExtraPath('inputImagesExp.xmd')
        if getSize(expImgMd) == 0:
            return

        level = self.level
        refImgMd = self.imgsRef
        self.classifyWholeSetStep(refImgMd, expImgMd, level)

        self.generateOutputMD(level)
        self.generateOutputClasses(level)

        self.levelUp(expImgMd)


    def _generateInputMd(self):
        #metadataItem = md.MetaData(self.imgsExp)
        #mdSize = md.getSize(self.imgsExp)
        #print("Len metadataItem", mdSize)
        metadataInput = md.MetaData()
        print("Len self.listOfParticles", len(self.listOfParticles), self.listToProc[0])

        count=0
        print("self.listToProc",self.listToProc)
        #print("self.listOfParticles", self.listOfParticles)
        #rows = iterRows(metadataItem)
        #for item in rows:
        for i in range(0,len(self.listOfParticles)): #AJ new:
            part = self.listOfParticles[0]
            imgRow = XmippMdRow() #AJ new:
            particleToRow(part, imgRow) #AJ new:
            objId = part.getObjId() #AJ new:
            print("objId", objId)
            #objId = item.getValue(md.MDL_ITEM_ID)
            #objIdFunc = item.getObjId()
            #print("objId", objIdFunc, objId)
            #if objId == self.listOfParticles[0]: #AJ new: remove this line
            #AJ new: copy the loadInput code lines
            imgRow.writeToMd(metadataInput, metadataInput.addObject())
            self.htAlreadyProcessed.pushItem(objId)
            count += 1
            self.listOfParticles.pop(0)
            if count==self.listToProc[0]:
                break

        self.numberProcessed += self.listToProc[0]
        self.totalInList -= self.listToProc[0]
        self.listToProc.pop(0)
        print("En _generateInputMd, Last. ", objId)
        metadataInput.write(self._getExtraPath('inputImagesExp.xmd')) #, md.MD_APPEND)

    def generateOutputClasses(self, level):

        if not exists(self._getExtraPath('last_classes.xmd')):

            finalMetadata = self._getExtraPath('final_classes.xmd')
            newMetadata = self._getExtraPath(join('level%04d' % level,
                                                  'general_level%04d' % level +
                                                  '_classes.xmd'))
            mdNew = md.MetaData('classes@' + newMetadata)
            if self.useAsRef == REF_CLASSES:
                mdRef = md.MetaData('classes@' + self.imgsRef)
            else:
                mdRef = md.MetaData(self.imgsRef)
            numList = []
            mdAll = md.MetaData()
            for item, itemRef in zip(mdNew, mdRef):
                imageAux = mdRef.getValue(md.MDL_IMAGE, itemRef)
                numAux = mdNew.getValue(md.MDL_REF, item)
                particles_total = mdNew.getValue(md.MDL_CLASS_COUNT, item)
                numList.append(numAux)
                row = md.Row()
                row.setValue(md.MDL_REF, numAux)
                row.setValue(md.MDL_IMAGE, imageAux)
                row.setValue(md.MDL_CLASS_COUNT, particles_total)
                row.addToMd(mdAll)
            mdAll.write('classes@' + finalMetadata, md.MD_APPEND)

            total = self.num_classes
            for i in range(total):
                self._params = {'newMd': 'class%06d_images@' % numList[i] +
                                         newMetadata,
                                'outMd': 'class%06d_images@' % numList[i] +
                                         finalMetadata}
                args = ('-i %(newMd)s -o %(outMd)s --mode append')
                self.runJob("xmipp_metadata_utilities",
                            args % self._params, numberOfMpi=1)

            copy(self._getExtraPath('final_classes.xmd'),
                 self._getExtraPath('last_classes.xmd'))

        else:
            finalMetadata = self._getExtraPath('final_classes.xmd')
            lastMetadata = self._getExtraPath('last_classes.xmd')
            newMetadata = self._getExtraPath(join('level%04d' % level,
                                                  'general_level%04d' % level +
                                                  '_classes.xmd'))

            mdNew = md.MetaData('classes@'+newMetadata)
            mdLast = md.MetaData('classes@'+lastMetadata)
            if self.useAsRef==REF_CLASSES:
                mdRef = md.MetaData('classes@'+self.imgsRef)
                total = getSize('classes@' + self.imgsRef)
            else:
                mdRef = md.MetaData(self.imgsRef)
                total = getSize(self.imgsRef)
            numList=[]
            mdAll = md.MetaData()
            for item, itemNew, itemRef in zip(mdNew, mdLast, mdRef):
                particles_total = mdLast.getValue(md.MDL_CLASS_COUNT, item) + \
                                  mdNew.getValue(md.MDL_CLASS_COUNT, itemNew)
                imageAux = mdRef.getValue(md.MDL_IMAGE, itemRef)
                numAux = mdNew.getValue(md.MDL_REF, item)
                numList.append(numAux)
                row = md.Row()
                row.setValue(md.MDL_REF, numAux)
                row.setValue(md.MDL_IMAGE, imageAux)
                row.setValue(md.MDL_CLASS_COUNT, particles_total)
                row.addToMd(mdAll)
            mdAll.write('classes@' + finalMetadata, md.MD_APPEND)

            for i in range(total):
                self._params = {'lastMd': 'class%06d_images@' % numList[i] +
                                          lastMetadata,
                                'newMd': 'class%06d_images@' % numList[i] +
                                         newMetadata,
                                'outMd': 'class%06d_images@' % numList[i] +
                                         finalMetadata}
                args = ('-i %(lastMd)s --set union_all %(newMd)s '
                        '-o %(outMd)s --mode append')
                self.runJob("xmipp_metadata_utilities",
                            args % self._params, numberOfMpi=1)

            copy(self._getExtraPath('final_classes.xmd'),
                 self._getExtraPath('last_classes.xmd'))


    def generateOutputMD(self, level):

        if not exists(self._getExtraPath('last_images.xmd')):
            copy(self._getExtraPath(join('level%04d' % level,
                                         'general_images_level%04d' % level +
                                         '.xmd')),
                 self._getExtraPath('last_images.xmd'))
            return

        self._params = {'lastMd': self._getExtraPath('last_images.xmd'),
                        'newMd': self._getExtraPath(join('level%04d' % level,
                                                         'general_images_level%04d'
                                                         % level + '.xmd')),
                        'outMd': self._getExtraPath('final_images.xmd')}
        args = ('-i %(lastMd)s --set union %(newMd)s -o %(outMd)s')
        self.runJob("xmipp_metadata_utilities",
                    args % self._params, numberOfMpi=1)

        copy(self._getExtraPath('final_images.xmd'),
             self._getExtraPath('last_images.xmd'))



    def classifyWholeSetStep(self, refImgMd, expImgMd, level):

        i=0
        while i <self.numberOfClassifyIterations:
            self.iterationStep(refImgMd, expImgMd, level)
            refImgMd = self._getExtraPath(join('level%04d' % level,
                                               'general_level%04d' % level +
                                               '_classes.xmd'))

            if i+1 < self.numberOfClassifyIterations:
                finalMetadata = self._getExtraPath(
                    join('level%04d' % level, 'general_level%04d' % level +
                         '_classes.xmd'))
                if exists(self._getExtraPath('final_classes.xmd')):
                    lastMetadata = self._getExtraPath('final_classes.xmd')
                else:
                    if level<2:
                        lastMetadata = self._getExtraPath('last_classes.xmd')
                    else:
                        lastMetadata = self._getExtraPath(
                            join('level%04d' % (level-1), 'general_level%04d'
                                 % (level-1) + '_classes.xmd'))
                newMetadata = self._getExtraPath(
                    join('level%04d' % level, 'general_level%04d' % level +
                         '_classes.xmd'))
                self.averageClasses(finalMetadata, lastMetadata, newMetadata, True)
            i += 1


    def iterationStep (self, refSet, imgsExp, level):

        if not exists(join(self._getExtraPath(), 'level%04d' % level)):
            mkdir(join(self._getExtraPath(), 'level%04d' % level))

        # Calling program xmipp_cuda_correlation
        filename = 'general_level%04d' % level + '_classes.xmd'
        self._params = {'imgsRef': refSet,
                        'imgsExp': imgsExp,
                        'outputFile': 'general_images_level%04d' % level + '.xmd',
                        'tmpDir': join(self._getExtraPath(),'level%04d' % level),
                        'keepBest': self.keepBest.get(),
                        'maxshift': self.maximumShift.get(),
                        'outputClassesFile': filename,
                        }

        args = '-i_ref %(imgsRef)s -i_exp %(imgsExp)s -o %(outputFile)s '\
               '--odir %(tmpDir)s --keep_best %(keepBest)d '\
               '--maxShift %(maxshift)d --simplifiedMd ' \
               '--classify %(outputClassesFile)s'
        self.runJob("xmipp_cuda_correlation", args % self._params, numberOfMpi=1)



    def createOutputStep(self):
        pass


    # --------------------------- UTILS functions -------------------------------

    def divideInBlocks(self):
        numOfImages = len(self.listOfParticles) - self.totalInList
        numOfBlk, lastBlk = divmod(numOfImages,self.blockToProc)
        self.listToProc += [self.blockToProc]*int(numOfBlk)
        if lastBlk != 0:
            self.listToProc.append(lastBlk)
            stepsToAppend = int(numOfBlk + 1)
        else:
            stepsToAppend = int(numOfBlk)
        self.totalInList += self.blockToProc*numOfBlk + lastBlk

        return stepsToAppend


    def _loadInputList(self):
        """ Load the input set of ctfs and create a list. """
        initial_time = time.time()
        particlesSet = self._loadInputParticleSet()
        self.isStreamClosed = particlesSet.getStreamState()
        self.listOfParticles = []
        for m in particlesSet:
            #imgRow = XmippMdRow()
            #particleToRow(m, imgRow)
            #idx = imgRow.getValue(md.MDL_ITEM_ID)
            idx = m.getObjId()
            if not self.htAlreadyProcessed.isItemPresent(idx):
                #self.listOfParticles.append(idx)
                newPart = m.clone()  # AJ new to make the list of particles instead of list of ids
                self.listOfParticles.append(newPart)
                #To create a md with just the new particles
                #imgRow.writeToMd(imagesMd, imagesMd.addObject())

        #imagesMd.write(self.imgsExp)

        particlesSet.close()
        self.debug("Closed db.")
        final_time = time.time()
        print("_loadInputList exec time",final_time-initial_time)


    def _loadInputParticleSet(self):
        initial_time = time.time()
        particlesFile = self.inputParticles.get().getFileName()
        self.debug("Loading input db: %s" % particlesFile)
        particlesSet = SetOfParticles(filename=particlesFile)
        particlesSet.loadAllProperties()
        final_time = time.time()
        print("_loadInputParticleSet exec time", final_time - initial_time)
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
        with open(AcceptedFile, 'w') as f: #AJ antes: 'a'
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

        #if exists(setFile):
        #    pwutils.path.cleanPath(setFile)

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


    def _updateOutputSetOfClasses(self, outputName, outputSet, state=Set.STREAM_OPEN):

        outputSet.setStreamState(state)
        if self.hasAttribute(outputName):
            print("En el if")
            outputSet.enableAppend()
            self._fillClassesFromLevel(outputSet, False)
            outputSet.write()  # Write to commit changes
            outputAttr = getattr(self, outputName)
            # Copy the properties to the object contained in the protocol
            outputAttr.copy(outputSet, copyId=False)
            # Persist changes
            self._store(outputAttr)
        else:
            print("En el else")
            self._fillClassesFromLevel(outputSet, True)
            self._defineOutputs(**{outputName: outputSet})
            self._defineSourceRelation(self.inputParticles, outputSet)
            self._store(outputSet)

        # Close set database to avoid locking it
        outputSet.close()

    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(md.MDL_REF))
        item.setTransform(rowToAlignment(row, ALIGN_2D))

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn, _, _ = self._classesInfo[classId]
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
            classCount = row.getValue(md.MDL_CLASS_COUNT)
            self._classesInfo[index] = (index, fn, classCount, row.clone())


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



    def levelUp(self, expImgMd):
        if getSize(expImgMd) == 0:
            return
        self.level += 1


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

