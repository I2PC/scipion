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

from pyworkflow.em.protocol import ProtAlign2D
import pyworkflow.em.metadata as md
import pyworkflow.protocol.params as params
from pyworkflow.em.metadata.utils import iterRows, getSize
from xmipp import MD_APPEND
from pyworkflow.em.packages.xmipp3.convert import rowToAlignment, \
    xmippToLocation
from convert import writeSetOfParticles, writeSetOfClasses2D
from shutil import copy
from os.path import join, exists
from os import mkdir, remove, listdir
import pyworkflow.protocol.constants as const
from pyworkflow.em import SetOfClasses2D, ALIGN_2D, ALIGN_NONE


class XmippProtGpuCrrCL2D(ProtAlign2D):
    """ 2D alignment using Xmipp GPU Correlation algorithm. """
    _label = 'gl2d'


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
                      pointerClass='SetOfParticles, SetOfAverages, '
                                   'SetOfClasses2D',
                      allowsNull=True,
                      label="Reference images",
                      help='Set of images that will serve as class reference')
        form.addParam('numberOfClasses', params.IntParam, default=32,
                      label='Number of classes:', important=True,
                      help='Number of classes (or references) to be generated')
        form.addParam('maxShift', params.IntParam, default=10,
                      label='Maximum shift (%):',
                      help='Maximum shift allowed during the alignment as '
                           'percentage of the input set size',
                      expertLevel=const.LEVEL_ADVANCED)
        form.addParam('keepBest', params.IntParam, default=1,
                      label='Number of best images:',
                      help='Number of classes to assign every input image '
                           'during the alignment',
                      expertLevel=const.LEVEL_ADVANCED)
        form.addParam('numberOfSplitIterations', params.IntParam, default=3,
                      label='Number of iterations in split stage:',
                      help='Maximum number of iterations in split stage',
                      expertLevel=const.LEVEL_ADVANCED)
        form.addParam('numberOfClassifyIterations', params.IntParam, default=5,
                      label='Number of iterations in classify stage:',
                      help='Maximum number of iterations when the  '
                           'classification of the whole image set is carried '
                           'out',
                      expertLevel=const.LEVEL_ADVANCED,)
        form.addParam('useAttraction', params.BooleanParam, default=True,
                      label='Allow attraction ?',
                      help='If you set to *Yes*, you allow to generate classes '
                           'with low number of images associated.\n'
                           'If *No*, all the generated classes will be '
                           'balanced',
                      expertLevel=const.LEVEL_ADVANCED,)
        form.addParallelSection(threads=0, mpi=0)


    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        """Mainly prepare the command line for call cuda corrrelation program"""

        # Convert input images if necessary
        if self.useReferenceImages:
            self.refSet = self._getExtraPath('imagesRef.xmd')
        else:
            self.refSet = None

        xOrig = self.inputParticles.get().getXDim()
        self.maximumShift = int(self.maxShift.get()*xOrig/100)
        self.p = 0.2
        self.percentStopClassify = 5
        self.listContinueClass=[]
        self.iterReturnSplit = 0
        self.iterReturnClass = 0
        self.depthSplit = 0
        self.depth = 0
        self.imgsExp = self._getExtraPath('imagesExp.xmd')

        self._insertFunctionStep('convertSetStep')
        self._insertFunctionStep('classifyStep')
        self._insertFunctionStep('createOutputStep')



    # --------------------------- STEPS functions --------------------------
    def convertSetStep(self):
        writeSetOfParticles(self.inputParticles.get(), self.imgsExp,
                            alignType=ALIGN_NONE)
        if self.useReferenceImages:
            if isinstance(self.referenceImages.get(), SetOfClasses2D):
                writeSetOfClasses2D(self.referenceImages.get(), self.refSet,
                                    writeParticles=False)
            else:
                writeSetOfParticles(self.referenceImages.get(), self.refSet,
                                alignType=ALIGN_NONE)


    def classifyStep(self):

        listNumImgs = []
        listNameImgs = []
        change = False

        level=0
        while len(listNumImgs) is not self.numberOfClasses.get() \
                or change is True:

            if self.useReferenceImages and self.numberOfClasses.get() \
                    == getSize(self.refSet):
                if not exists(join(self._getExtraPath(), 'level%03d' % level)):
                    mkdir(join(self._getExtraPath(), 'level%03d' % level))
                copy(self.refSet, self._getExtraPath(join('level%03d' % level,
                                                'intermediate_classes.xmd')))
            else:
                self.splitStep(level)

                #######################################
                if not self.useAttraction:
                    self.attractionSplitStep(level)
                #######################################

                self.generateMetadata(listNameImgs, listNumImgs, level)

            self.ref = self._getExtraPath(join('level%03d' % level,
                                               'intermediate_classes.xmd'))
            lengthMd = getSize(self.ref)
            if lengthMd==self.numberOfClasses.get():
                self.classifyWholeSetStep(level)
                #######################################################
                if not self.useAttraction:
                    change, listNumImgs, listNameImgs = \
                        self.attractionGeneralStep(level)
                    if change:
                        level = level + 1
                        continue
                #######################################################
                copy(self._getExtraPath(join('level%03d' % level,
                                             'general_level%03d' % level
                                             + '_classes.xmd')),
                     self._getExtraPath('last_classes.xmd'), )
                if exists(self._getExtraPath(join('level%03d' % level,
                                                  'general_images_level%03d'
                                                  % level + '.xmd'))):
                    copy(self._getExtraPath(join('level%03d' % level,
                                                 'general_images_level%03d'
                                                 % level + '.xmd')),
                         self._getExtraPath('last_images.xmd'), )
                return

            listNumImgs, listNameImgs = self.checkOutput(level)
            self.cleaningPath(level)
            level = level + 1


    def createOutputStep(self):

        inputParticles = self.inputParticles.get()
        classes2DSet = self._createSetOfClasses2D(inputParticles)
        self._fillClassesFromLevel(classes2DSet)
        result = {'outputClasses': classes2DSet}
        self._defineOutputs(**result)
        self._defineSourceRelation(self.inputParticles, classes2DSet)


    # --------------------------- Other functions --------------------------

    def splitStep(self, level):

        expSet = self.imgsExp
        if level > 0:
            expSet = self._getExtraPath(join('level%03d' % (level-1),
                                             'images_level%03d' % (level-1)
                                             + '_major.xmd'))
            boolReferenceImages = False
        else:
            boolReferenceImages = self.useReferenceImages

        i=0
        while i <self.numberOfSplitIterations:
            self.iterationStep(self.refSet, expSet, i, boolReferenceImages,
                               level, True, False)
            self.refSet = self._getExtraPath(join('level%03d' % level,
                                                  'level%03d' % level +
                                                  '_classes.xmd'))
            i+=1
            if not self.useAttraction:
                if self.fastCheckingAtt(level, True):
                    self.iterReturnSplit=i
                    return


    def classifyWholeSetStep(self, level):

        i=self.iterReturnClass
        if i == self.numberOfClassifyIterations:
            copy(self._getExtraPath(join('level%03d' % level,
                                         'intermediate_classes.xmd')),
                 self._getExtraPath(join('level%03d' % level,
                                         'general_level%03d_classes' % level
                                         +  '.xmd')))

        while i <self.numberOfClassifyIterations:
            self.iterationStep(self.ref, self.imgsExp, i, False, level, False,
                               False)
            self.ref = self._getExtraPath(join('level%03d' % level,
                                               'general_level%03d' % level +
                                               '_classes.xmd'))
            i+=1

            if self.checkContinueClassification(level, i-1):
                return

            #check attraction
            if not self.useAttraction:
                if self.fastCheckingAtt(level, False):
                    self.iterReturnClass=i
                    return


    #check attraction
    def fastCheckingAtt (self, level, flag_split):
        listAuxNum = []
        if flag_split:
            metadata = md.MetaData(self._getExtraPath(join('level%03d' % level,
                                   'level%03d' % level + '_classes.xmd')))
        else:
            metadata = md.MetaData(self._getExtraPath(join('level%03d' % level,
                                'general_level%03d' % level + '_classes.xmd')))
        for item in metadata:
            numImgs = metadata.getValue(md.MDL_CLASS_COUNT, item)
            listAuxNum.append(numImgs)
        total = sum(listAuxNum)
        th = (self.p * total / len(listAuxNum))
        aux = [i for i in listAuxNum if i<th]
        if len(aux)>0:
            return True
        else:
            return False


    def checkContinueClassification(self, level, iter):

        diff=0
        i=0
        metadata = md.MetaData(self._getExtraPath(join('level%03d' % level,
                                'general_images_level%03d' % level + '.xmd')))

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
        if num<self.percentStopClassify and iter>0:
            return True
        else:
            return False


    def iterationStep (self, refSet, imgsExp, iter, useReferenceImages,
                       level, flag_split, flag_attraction):

        if not exists(join(self._getExtraPath(), 'level%03d' % level)):
            mkdir(join(self._getExtraPath(), 'level%03d' % level))

        if not useReferenceImages and iter==0 and flag_split==True:

            # First step: divide the metadata input file to generate
            # a couple of references
            if level==0:
                if not flag_attraction:
                    outDirName = imgsExp[0:imgsExp.find('extra')+6] +  \
                                 'level%03d' % level +  imgsExp[imgsExp.find(
                        'extra')+5:-4]
                else:
                    outDirName = imgsExp[0:imgsExp.find('extra') + 6] +  \
                                 'level%03d' % level +  imgsExp[imgsExp.find(
                        'level%03d' % level) + 8:-4]
            else:
                if not flag_attraction:
                    outDirName = imgsExp[0:imgsExp.find('extra') + 6] +  \
                                 'level%03d' % level +  imgsExp[imgsExp.find(
                        'level%03d' % (level-1)) + 8:-4]
                else:
                    outDirName = imgsExp[0:imgsExp.find('extra') + 6] +  \
                                 'level%03d' % level +  imgsExp[imgsExp.find(
                        'level%03d' % level) + 8:-4]
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

            # Third step: generate a single metadata with the two previous
            # averages
            refSet = self._getExtraPath(join('level%03d' % level,'refSet.xmd'))
            self._params = {'avg1': avg1 + 'average.xmp',
                            'avg2': avg2 + 'average.xmp',
                            'outputMd': refSet}
            args = ('-i %(avg1)s --set union %(avg2)s -o %(outputMd)s')
            self.runJob("xmipp_metadata_utilities", args % self._params,
                        numberOfMpi=1)

        # Fourth step: calling program xmipp_cuda_correlation
        if flag_split:
            filename = 'level%03d' % level+'_classes.xmd'
            self._params = {'imgsRef': refSet,
                            'imgsExp': imgsExp,
                            'outputFile': 'images_level%03d' % level+'.xmd',
                            'tmpDir': join(self._getExtraPath(),'level%03d'
                                           % level),
                            'keepBest': self.keepBest.get(),
                            'maxshift': self.maximumShift,
                            'outputClassesFile': filename,
                            }
        else:
            filename = 'general_level%03d' % level + '_classes.xmd'
            self._params = {'imgsRef': refSet,
                            'imgsExp': imgsExp,
                            'outputFile': 'general_images_level%03d' % level
                                          + '.xmd',
                            'tmpDir': join(self._getExtraPath(),'level%03d'
                                           % level),
                            'keepBest': self.keepBest.get(),
                            'maxshift': self.maximumShift,
                            'outputClassesFile': filename,
                            }
        Nrefs = getSize(refSet)
        if Nrefs>2:
            args = '-i_ref %(imgsRef)s -i_exp %(imgsExp)s -o %(outputFile)s '\
                   '--odir %(tmpDir)s --keep_best %(keepBest)d '\
                   '--maxShift %(maxshift)d --classify %(outputClassesFile)s '\
                   '--simplifiedMd'
            self.runJob("xmipp_cuda_correlation", args % self._params,
                        numberOfMpi=1)
        else:
            flag_error=True
            while flag_error:
                try:
                    flag_error=False
                    args='-i %(imgsExp)s --ref0 %(imgsRef)s --nref %(Nrefs)d '\
                         '--iter 1 --distance correlation --classicalMultiref '\
                         '--maxShift %(maxshift)d --odir %(cl2dDir)s'
                    self._params['Nrefs']=Nrefs
                    self._params['cl2dDir'] = self._getExtraPath(
                                              join('level%03d' % level))
                    self.runJob("mpirun -np 4 -bynode xmipp_mpi_classify_CL2D",
                                args % self._params)
                except Exception as ex:
                    flag_error = True

            if flag_split:
                copy(self._getExtraPath(join('level%03d' % level,
                                             "level_00","class_classes.xmd")),
                     self._getExtraPath(join('level%03d' % level,
                                             'level%03d' % level
                                             + '_classes.xmd')))
                copy(self._getExtraPath(join('level%03d' % level,"images.xmd")),
                     self._getExtraPath(join('level%03d' % level,
                                             'images_level%03d' % level +
                                             '.xmd')))
            else:
                copy(self._getExtraPath(join('level%03d' % level,
                                             "level_00", "class_classes.xmd")),
                     self._getExtraPath(join('level%03d' % level,
                                             'general_level%03d' % level +
                                             '_classes.xmd')))
                copy(self._getExtraPath(join('level%03d' % level,
                                             "images.xmd")),
                     self._getExtraPath(join('level%03d' % level,
                                             'general_images_level%03d' %
                                             level + '.xmd')))


    def attractionSplitStep(self, level):

        change, labelMaxClass, labelMinClass, mdToReduce, mdToCheck, \
        listNumImgs, listNameImgs = self.checkAttraction(level, True)

        while change:

            #Every time we need to make a change (for unbalanced classes) we
            # update with +1 this value
            self.depthSplit += 1

            self.imgsExp = self._getExtraPath(join('level%03d' % level,
                                                   'images_level%03d' % level
                                                   + '_NoAtt.xmd'))
            self._params = {'input': 'class%06d_images' % (labelMaxClass)  +
                                     '@' + mdToCheck,
                            'outputMd': self.imgsExp}
            args = ('-i %(input)s -o %(outputMd)s')
            self.runJob("xmipp_metadata_utilities", args % self._params,
                        numberOfMpi=1)

            i = 0
            while i < self.numberOfSplitIterations:
                self.iterationStep(self.refSet, self.imgsExp, i, False,
                                   level, True, True)
                self.refSet = self._getExtraPath(join('level%03d' % level,
                                                      'level%03d' % level +
                                                      '_classes.xmd'))
                i+=1
                if self.fastCheckingAtt(level, True):
                    break

            #Recursive calling of this function
            self.attractionSplitStep(level)

            #When we detect the classes are balanced, we update with -1 this
            # value, because we finish one recursive calling of the function
            if self.depthSplit>1:
                self.depthSplit-=1
                return
            if self.depthSplit==1:
                self.depthSplit=0

            if (level - 1) >= 0:
                self.imgsExp = self._getExtraPath(join('level%03d' % (
                        level-1), 'images_level%03d'%(level-1) + '_major.xmd'))
            else:
                self.imgsExp = self._getExtraPath('imagesExp.xmd')

            i = self.iterReturnSplit
            if i==self.numberOfSplitIterations:
                i-=1
            while i < self.numberOfSplitIterations:
                self.iterationStep(self.refSet, self.imgsExp, i, True, level,
                                   True, False)
                self.refSet = self._getExtraPath(join('level%03d' % level,
                                                      'level%03d' % level +
                                                      '_classes.xmd'))
                i+=1
                if self.fastCheckingAtt(level, True):
                    break

            change, labelMaxClass, __, mdToReduce, mdToCheck, \
            __, __ = self.checkAttraction(level, True)



    def attractionGeneralStep(self, level):

        change, labelMaxClass, labelMinClass, mdToReduce, mdToCheck, \
        listNumImgs, listNameImgs = self.checkAttraction(level, False)

        if change:
            self.depth+=1

        if change:
            self._params = {'input': 'class%06d_images' % (labelMaxClass)  +
                                     '@' + mdToCheck,
                            'outputMd': self._getExtraPath(join('level%03d'
                                        % level, 'images_level%03d' % level
                                                                + '_major.xmd'))
                            }
            args = ('-i %(input)s -o %(outputMd)s')
            self.runJob("xmipp_metadata_utilities", args % self._params,
                        numberOfMpi=1)

        return change, listNumImgs, listNameImgs



    def checkAttraction(self, level, flag_split):

        if flag_split:
            mdToCheck = self._getExtraPath(join('level%03d' % level,
                                                'level%03d' % level
                                                + '_classes.xmd'))
            mdToReduce = self._getExtraPath(join('level%03d' % level,
                                                 'images_level%03d' % level
                                                 + '.xmd'))
        else:
            mdToCheck = self._getExtraPath(join('level%03d' % level,
                                                'general_level%03d' % level
                                                + '_classes.xmd'))
            mdToReduce = self._getExtraPath(join('level%03d' % level,
                                                 'general_images_level%03d'
                                                 % level + '.xmd'))

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
        total = sum(listAuxNum)

        th = (self.p*total/len(listAuxNum))
        labelMinClass=[]
        for i in range(len(listAuxNum)):
            if listAuxNum[i]<th:
                labelMinClass.append(listAuxRef[i])
                listAuxNum[i] = -1

        labelMaxClass = listAuxRef[listAuxNum.index(max(listAuxNum))]
        listAuxNum[labelMaxClass-1] = -1

        if len(labelMinClass)>0:
            change = True
        else:
            change = False

        return change, labelMaxClass, labelMinClass, mdToReduce, mdToCheck, \
               listAuxNum, listAuxName


    def checkOutput(self, level):

        listAuxString = []
        listAuxNum = []
        listAuxRefs = []

        outSet = self._getExtraPath(join('level%03d' % level,
                                         'intermediate_classes.xmd'))
        metadataItem = md.MetaData(outSet)
        for item in metadataItem:
            nameImg = metadataItem.getValue(md.MDL_IMAGE, item)
            listAuxString.append(nameImg)
            numImgs = metadataItem.getValue(md.MDL_CLASS_COUNT, item)
            listAuxNum.append(numImgs)
            refImg = metadataItem.getValue(md.MDL_REF, item)
            listAuxRefs.append(refImg)

        maxValue = max(listAuxNum)
        maxPos = listAuxNum.index(maxValue)

        listAuxNum[maxPos] = -1
        bestRef = listAuxRefs[maxPos]

        self._params = {'input': 'class%06d_images' % (bestRef) + '@' + outSet,
                        'outputMd': self._getExtraPath(join('level%03d'  %
                                    level,'images_level%03d' % level
                                                            + '_major.xmd'))
                        }
        args = ('-i %(input)s -o %(outputMd)s')
        self.runJob("xmipp_metadata_utilities", args % self._params,
                    numberOfMpi=1)

        return listAuxNum, listAuxString


    def generateMetadata(self, listNameImgs, listNumImgs, level):

        # Renumerate unchanged classes
        listNewNumImages = [-1] * len(listNumImgs)
        count = 1
        for i in range(len(listNumImgs)):
            if listNumImgs[i] is not -1:
                listNewNumImages[i] = count
                count += 1

        # Construct the new classes with the renumerated old classes
        mdNewClasses = md.MetaData()
        for i in range(len(listNumImgs)):
            if listNumImgs[i] is not -1:
                name = listNameImgs[i]
                fn = name[name.find('@') + 1:-4] + '.xmd'
                numRef = int(name[0:6])

                mdClass = md.MetaData("classes@" + fn)
                for row in iterRows(mdClass):
                    if mdClass.getValue(md.MDL_REF, row.getObjId()) == numRef:
                        row.setValue(md.MDL_REF, listNewNumImages[i])
                        row.addToMd(mdNewClasses)

        # Add the two new classes to the list of renumerated classes
        outSet = self._getExtraPath(join('level%03d' % level, 'level%03d' %
                                         level + '_classes.xmd'))
        mdClass = md.MetaData("classes@" + outSet)
        rows = iterRows(mdClass)
        for row in rows:
            row.setValue(md.MDL_REF, count)
            row.addToMd(mdNewClasses)
            count = count + 1
        mdNewClasses.write('classes@' + self._getExtraPath(join('level%03d'
                                        % level,'intermediate_classes.xmd')),
                           MD_APPEND)

        # Generate the intermediate images and the blocks of the intermediate
        # classes for the unchanged classes
        for i in range(len(listNumImgs)):
            if listNumImgs[i] is not -1:
                name = listNameImgs[i]
                fn = name[name.find('@') + 1:-4] + '.xmd'
                numRef = int(name[0:6])

                # Read the list of images in this class
                mdImgsInClass = md.MetaData(
                    'class%06d_images@%s' % (numRef, fn))
                mdImgsInClass.fillConstant(md.MDL_REF, listNewNumImages[i])
                mdImgsInClass.write('class%06d' % (listNewNumImages[i])
                                    + '_images@' + self._getExtraPath(join(
                                    'level%03d' % level,
                                    'intermediate_classes.xmd')),
                                    MD_APPEND)

        # Add the two new classes
        if len(listNumImgs) == 0:
            count = 1
        else:
            count = len(listNumImgs)
        for newRef in range(getSize(outSet)):
            mdImgsInClass = md.MetaData('class%06d_images@%s' % (newRef + 1,
                                                                 outSet))
            mdImgsInClass.fillConstant(md.MDL_REF, count)
            mdImgsInClass.write('class%06d' % (count) + '_images@'
                                + self._getExtraPath(join('level%03d' %
                                level,'intermediate_classes.xmd')),
                                MD_APPEND)
            count = count + 1



    def cleaningPath(self, level):
        if exists(self._getExtraPath(join('level%03d' % level,'refSet.xmd'))):
            remove(self._getExtraPath(join('level%03d' % level,'refSet.xmd')))
        if exists(self._getExtraPath(join('level%03d' % level,"images.xmd"))):
            remove(self._getExtraPath(join('level%03d' % level,"images.xmd")))
            level_1=level-1
            if level>0 and exists(self._getExtraPath(join('level%03d' %
                            level_1,"images_level%03d_major.xmd" %level_1))):
                remove(self._getExtraPath(join('level%03d' % level_1,
                                               "images_level%03d_major.xmd"
                                               %level_1)))
        for f in listdir(join(join(self._getExtraPath(),'level%03d'%level))):
            if not f.find('average')==-1 \
                    or not f.find('stddev')==-1 \
                    or not f.find('NoAtt')==-1:
                remove(join(join(self._getExtraPath(),'level%03d'%level,f)))
            if level>0 and not f.find('major0') == -1:
                remove(join(join(self._getExtraPath(), 'level%03d' % level, f)))
            if level==0 and not f.find('imagesExp0') == -1:
                remove(join(join(self._getExtraPath(), 'level%03d' % level, f)))



    # --------------------------- UTILS functions ------------------------------

    def _fillClassesFromLevel(self, clsSet):
        """ Create the SetOfClasses2D from a given iteration. """
        myFileClasses = self._getExtraPath('last_classes.xmd')
        myFileImages = self._getExtraPath('last_images.xmd')
        self._loadClassesInfo(myFileClasses)
        xmpMd = myFileImages
        iterator = md.SetMdIterator(xmpMd, sortByLabel=md.MDL_ITEM_ID,
                                    updateItemCallback=self._updateParticle,
                                    skipDisabled=True)
        clsSet.classifyItems(updateItemCallback=iterator.updateItem,
                             updateClassCallback=self._updateClass)

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


    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        if self.useReferenceImages:
            if self.numberOfClasses<self.referenceImages.get().getSize():
                errors.append('The number of classes must be equal or greater'
                              ' than the number of references')
            if self.referenceImages.hasValue():
                [x1, y1, z1] = self.referenceImages.get().getDimensions()
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

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputClasses'):
            methods.append("Output alignment not ready yet.")
        else:
            if self.useReferenceImages:
                methods.append("We aligned images %s with respect to the "
                               "reference image set %s using Xmipp GPU  "
                               "correlation method"
                    % (self.getObjectTag('inputParticles'),
                       self.getObjectTag('referenceImages')))
            else:
                methods.append(
                    "We aligned images %s with no references using Xmipp GPU "
                    "correlation method" % self.getObjectTag('inputParticles'))
            methods.append(" and produced %s images."
                           % self.getObjectTag('outputClasses'))
        return methods

