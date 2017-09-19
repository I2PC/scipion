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

from pyworkflow.em import ALIGN_NONE, ALIGN_2D
from pyworkflow.em.protocol import ProtAlign2D
import pyworkflow.em.metadata as md
import pyworkflow.protocol.params as params
from pyworkflow.em.metadata.utils import iterRows, getSize
from xmipp import MD_APPEND
from pyworkflow.em.packages.xmipp3.convert import readSetOfClasses2D
from convert import writeSetOfParticles
from shutil import copy


class XmippProtGpuCrrCL2D(ProtAlign2D):
    """ Aligns a set of particles using the GPU Correlation algorithm. """
    _label = 'align with GPU Correlation'


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
        form.addParam('keepBest', params.IntParam, default=1,
                      label='Number of best images:',
                      help='Number of the best images to keep for every class')
        form.addParam('numberOfIterations', params.IntParam, default=2, #DEJAR A 10
                      label='Number of iterations:',
                      help='Maximum number of iterations')
        form.addParam('numberOfClasses', params.IntParam, default=5, #DEJAR A 64
                      label='Number of classes:',
                      help='Number of classes (or references) to be generated.')
        form.addParam('useAttraction', params.BooleanParam, default=True,
                      label='Allow attraction ?',
                      help='If you set to *Yes*, you allow to generate classes '
                           'with low number of images associated.\n'
                           'If *No*, all the generated classes will be '
                           'balanced.')


    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        """" Mainly prepare the command line for call cuda corrrelation program"""

        # Convert input images if necessary
        if self.useReferenceImages:
            self.refSet = self._getExtraPath('imagesRef.xmd')
            self._insertFunctionStep('convertSetStep', self.refSet, False)
        else:
            self.refSet = None

        self.levels = []
        self.numRefs = []
        self.cadena = []
        self.depthSplit = 0
        self.depth = 0
        self.imgsExp = self._getExtraPath('imagesExp.xmd')
        self.exp = self.imgsExp
        self._insertFunctionStep('convertSetStep', self.imgsExp, True)

        self._insertFunctionStep('classifyStep')

        self._insertFunctionStep('createOutputStep')



    # --------------------------- STEPS functions --------------------------
    def convertSetStep(self, imgs, exp):
        if exp:
            writeSetOfParticles(self.inputParticles.get(), imgs,
                                alignType=ALIGN_NONE)
        else:
            writeSetOfParticles(self.referenceImages.get(), imgs,
                                alignType=ALIGN_NONE)


    def generateMetadata(self, listNameImgs, listNumImgs, level):


        print('generateMetadata', level)
        print('listNameImgs en generateMetadata', listNameImgs)
        print('listNumImgs en generateMetadata', listNumImgs)

        # if(level>0):
        #
        #     testMd = self._getExtraPath('test_level00' + str(level-1) + '.xmd')
        #     self._params = {'input': testMd,
        #                     'numRef': listRefImgs[listNumImgs.index(max(listNumImgs))],
        #                     'outputMd': self._getExtraPath('test_level00' + str(level-1) + '_minor.xmd')}
        #     args = ('-i %(input)s --query select "ref==%(numRef)d" -o %(outputMd)s')
        #     self.runJob("xmipp_metadata_utilities", args % self._params)
        #
        #     testMd = self._getExtraPath('test_level00' + str(level-1) + '_minor.xmd')
        #     outMd = self._getExtraPath('intermediate_test.xmd')
        #     self._params = {'input': testMd,
        #                     'secondMd': self._getExtraPath('test_level00' + str(level) + '.xmd'),
        #                     'outputMd': outMd}
        #     args = ('-i %(input)s --set union %(secondMd)s -o %(outputMd)s')
        #     self.runJob("xmipp_metadata_utilities", args % self._params)
        # else:
        #     copy(self._getExtraPath('test_level00' + str(level) + '.xmd'), self._getExtraPath('intermediate_test.xmd'))


        mdAll = md.MetaData()
        count = 0
        for i in range(len(listNumImgs)):
            if listNumImgs[i] is not -1:
                name = listNameImgs[i]
                numRef = int(name[0:6])
                print(name)

                fn = name[name.find('@')+1:-4]+'.xmd'
                blocks = md.getBlocksInMetaDataFile(fn)

                for block in blocks:
                    if block.startswith('classes'):
                        mdClass = md.MetaData(block + "@" + fn)
                        rows = iterRows(mdClass)
                        for row in rows:
                            if mdClass.getValue(md.MDL_REF, row.getObjId()) == numRef:
                                row.setValue(md.MDL_REF, count)
                                row.addToMd(mdAll)
                                count = count + 1
                                break

        outSet = self._getExtraPath('level00' + str(level) + '_classes.xmd')
        blocks = md.getBlocksInMetaDataFile(outSet)
        for block in blocks:
            if block.startswith('classes'):
                mdClass = md.MetaData(block + "@" + outSet)
                rows = iterRows(mdClass)
                for row in rows:
                    row.setValue(md.MDL_REF, count)
                    row.addToMd(mdAll)
                    count = count + 1

        mdAll.write('classes@' + self._getExtraPath('intermediate_classes.xmd'), MD_APPEND)

        count = 0
        for i in range(len(listNumImgs)):
            if listNumImgs[i] is not -1:
                name = listNameImgs[i]
                numRef = int(name[0:6])

                fn = name[name.find('@')+1:-4]+'.xmd'
                mdAll2 = md.MetaData()
                blocks = md.getBlocksInMetaDataFile(fn)

                for block in blocks:
                    if block.startswith('class%06d' % (numRef)):
                        mdClass = md.MetaData(block + "@" + fn)
                        rows = iterRows(mdClass)
                        for row in rows:
                            row.addToMd(mdAll2)
                        mdAll2.write('class%06d' % (count) + '_images@' + self._getExtraPath('intermediate_classes.xmd'), MD_APPEND)
                        count = count + 1
                        break

        outSet = self._getExtraPath('level00' + str(level) + '_classes.xmd')
        blocks = md.getBlocksInMetaDataFile(outSet)
        for block in blocks:
            mdAll2 = md.MetaData()
            if block.startswith('class0'):
                mdClass = md.MetaData(block + "@" + outSet)
                rows = iterRows(mdClass)
                for row in rows:
                    row.addToMd(mdAll2)
                mdAll2.write('class%06d' % (count) + '_images@' + self._getExtraPath('intermediate_classes.xmd'), MD_APPEND)
                count = count + 1





    # def generateMetadata(self, levels, numRefs):
    #
    #     mdAll = md.MetaData()
    #     count=0
    #     for i in range(len(levels)):
    #         numRef = numRefs[i]
    #         level = levels[i]
    #
    #         fn = self._getExtraPath('level00'+str(level)+'_iter'+str(self.numberOfIterations.get()-1)+'_classes.xmd')
    #         blocks = md.getBlocksInMetaDataFile(fn)
    #
    #         for block in blocks:
    #             if block.startswith('classes'):
    #                 mdClass = md.MetaData(block + "@" + fn)
    #                 rows = iterRows(mdClass)
    #                 for row in rows:
    #                     if mdClass.getValue(md.MDL_REF, row.getObjId())==numRef+1:
    #                         row.setValue(md.MDL_REF, count)
    #                         row.addToMd(mdAll)
    #                         count=count+1
    #                         break
    #     mdAll.write('classes@' + self._getExtraPath('final_classes.xmd'), MD_APPEND)
    #
    #
    #     count=0
    #     for i in range(len(levels)):
    #         numRef = numRefs[i]
    #         level = levels[i]
    #
    #         fn = self._getExtraPath('level00' + str(level) + '_iter' + str(self.numberOfIterations.get() - 1) + '_classes.xmd')
    #         mdAll2 = md.MetaData()
    #         blocks = md.getBlocksInMetaDataFile(fn)
    #
    #         for block in blocks:
    #             if block.startswith('class%06d'%(numRef+1)):
    #                 mdClass = md.MetaData(block + "@" + fn)
    #                 rows = iterRows(mdClass)
    #                 for row in rows:
    #                     row.addToMd(mdAll2)
    #                 mdAll2.write('class%06d'%(count)+'_images@'+self._getExtraPath('final_classes.xmd'), MD_APPEND)
    #                 count=count+1
    #                 break


    def classifyStep(self):

        print('classifyStep')

        listNumImgs = []
        listNameImgs = []
        flag_repeat=False
        change = False

        level=0
        while len(listNumImgs) is not self.numberOfClasses.get() or change is True:

            self.splitStep(level, flag_repeat)

            #######################################
            if not self.useAttraction:
                self.attractionSplitStep(level)
            #######################################

            self.generateMetadata(listNameImgs, listNumImgs, level)
            self.ref = self._getExtraPath('intermediate_classes.xmd')

            lengthMd = getSize(self.ref)
            if lengthMd==self.numberOfClasses.get():
                self.classifyWholeSetStep(level)
                ##############################################################################
                if not self.useAttraction:
                    change, listNumImgs, listNameImgs = self.attractionGeneralStep(level)
                    if change:
                        level = level + 1
                        continue
                ##############################################################################
                copy(self._getExtraPath('general_level00' + str(level) + '_classes.xmd'), self._getExtraPath('last_classes.xmd'), )
                return
            else:
                copy(self._getExtraPath('intermediate_classes.xmd'), self._getExtraPath('general_level00'+str(level)+'_classes.xmd'))


            listNumImgs, listNameImgs = self.checkOutput(level)
            level = level + 1


    def splitStep(self, level, flag_repeat):

        print('splitStep')

        if level > 0:
            self.imgsExp = self._getExtraPath('test_level00' + str(level - 1) + '_major.xmd')
            boolReferenceImages = False
        else:
            boolReferenceImages = self.useReferenceImages

        i=0
        while i <self.numberOfIterations:
            self.iterationStep(self.refSet, self.imgsExp, i, boolReferenceImages, level, True)
            self.refSet = self._getExtraPath('level00' + str(level) + '_classes.xmd')
            i+=1

            # AJ comprobar si hay dos clases o no porque si no las hay, TENEMOS UN PROBLEMA
            length = getSize(self._getExtraPath('level00' + str(level) + '_classes.xmd'))
            if length == 1:
                print('PROBLEMAAAA')
                i = 0




    def classifyWholeSetStep(self, level):

        print('classifyWholeSetStep')

        boolReferenceImages = False

        for i in range(self.numberOfIterations):
            self.iterationStep(self.ref, self.exp, i, boolReferenceImages, level, False)
            self.ref = self._getExtraPath('general_level00' + str(level) + '_classes.xmd')



    def iterationStep (self, refSet, imgsExp, iter, useReferenceImages, level, flag_split):

        print('iterationStep')

        if not useReferenceImages and iter==0 and flag_split==True:
            # First step: divide the metadata input file to generate
            # a couple of references
            self._params = {'imgsExp': imgsExp}
            args = ('-i %(imgsExp)s -n 2')
            self.runJob("xmipp_metadata_split", args % self._params)

            # Second step: calculate the means of the previous metadata
            expSet1 = imgsExp[0:-4] + '000001.xmd'
            avg1 = imgsExp[0:-4] + '_000001'
            expSet2 = imgsExp[0:-4] + '000002.xmd'
            avg2 = imgsExp[0:-4] + '_000002'
            self._params = {'imgsSet': expSet1,
                            'outputAvg': avg1}
            args = ('-i %(imgsSet)s --save_image_stats %(outputAvg)s -v 0')
            self.runJob("xmipp_image_statistics", args % self._params)

            self._params = {'imgsSet': expSet2,
                            'outputAvg': avg2}
            args = ('-i %(imgsSet)s --save_image_stats %(outputAvg)s -v 0')
            self.runJob("xmipp_image_statistics", args % self._params)

            # Third step: generate a single metadata with the two previous averages
            refSet = self._getExtraPath('refSet.xmd')
            self._params = {'avg1': avg1 + 'average.xmp',
                            'avg2': avg2 + 'average.xmp',
                            'outputMd': refSet}
            args = ('-i %(avg1)s --set union %(avg2)s -o %(outputMd)s')
            self.runJob("xmipp_metadata_utilities", args % self._params)

        # Fourth step: calling program xmipp_cuda_correlation
        if flag_split:
            filename = 'level00'+str(level)+'_classes.xmd'
            self._params = {'imgsRef': refSet,
                            'imgsExp': imgsExp,
                            'outputFile': 'test_level00'+str(level)+'.xmd',
                            'tmpDir': self._getExtraPath(),
                            'keepBest': self.keepBest.get(),
                            'maxshift': self.maximumShift.get(),
                            'outputClassesFile': filename,
                            }
        else:
            filename = 'general_level00' + str(level) + '_classes.xmd'
            self._params = {'imgsRef': refSet,
                            'imgsExp': imgsExp,
                            'outputFile': 'general_test_level00' + str(level) + '.xmd',
                            'tmpDir': self._getExtraPath(),
                            'keepBest': self.keepBest.get(),
                            'maxshift': self.maximumShift.get(),
                            'outputClassesFile': filename,
                            }
        args = ('-i_ref %(imgsRef)s -i_exp %(imgsExp)s -o %(outputFile)s '
                '--odir %(tmpDir)s --keep_best %(keepBest)d '
                '--maxShift %(maxshift)d --classify %(outputClassesFile)s')
        self.runJob("xmipp_cuda_correlation", args % self._params)



    def attractionSplitStep(self, level):

        change, labelMaxClass, __, mdToReduce, __, __ = self.checkAttraction(level, True)
        if change:
            self.depthSplit+=1

        while change:
            print('CHANGEEEEEEEEEEEEEEEE')

            self._params = {'input': mdToReduce,
                            'numRef': labelMaxClass,
                            'outputMd': mdToReduce[0:-4]+'_NoAtt.xmd'}
            args = ('-i %(input)s --query select "ref==%(numRef)d" -o %(outputMd)s')
            self.runJob("xmipp_metadata_utilities", args % self._params)

            self.imgsExp = self._getExtraPath('test_level00' + str(level) + '_NoAtt.xmd')

            i = 0
            while i < self.numberOfIterations:
                self.iterationStep(self.refSet, self.imgsExp, i, False, level, True)
                self.refSet = self._getExtraPath('level00' + str(level) + '_classes.xmd')
                i+=1

                # AJ comprobar si hay dos clases o no porque si no las hay, TENEMOS UN PROBLEMA
                length = getSize(self._getExtraPath('level00' + str(level) + '_classes.xmd'))
                if length == 1:
                    print('PROBLEMAAAA')
                    i = 0

            self.attractionSplitStep(level)

            if self.depthSplit>1:
                self.depthSplit-=1
                return
            if self.depthSplit==1:
                self.depthSplit=0

            if (level - 1) >= 0:
                self.imgsExp = self._getExtraPath('test_level00' + str(level - 1) + '_major.xmd')
            else:
                self.imgsExp = self._getExtraPath('imagesExp.xmd')

            i = 0
            while i < self.numberOfIterations:
                self.iterationStep(self.refSet, self.imgsExp, i, True, level, True)
                self.refSet = self._getExtraPath('level00' + str(level) + '_classes.xmd')
                i+=1

                # AJ comprobar si hay dos clases o no porque si no las hay, TENEMOS UN PROBLEMA
                length = getSize(self._getExtraPath('level00' + str(level) + '_classes.xmd'))
                if length == 1:
                    print('PROBLEMAAAA')
                    i = 0

            change, labelMaxClass, __, mdToReduce, __, __ = self.checkAttraction(level, True)



    def attractionGeneralStep(self, level):

        change, labelMaxClass, labelMinClass, mdToReduce, listNumImgs, listNameImgs = self.checkAttraction(level, False)
        #change = True
        #for i in range(len(listNumImgs)):
        #    if listNumImgs[i]!=-1:
        #        listNumImgs[i]=-1
        #    if listNumImgs.count(-1)==3:
        #       break

        if change:
            self.depth+=1

        if change: #esto era un while
            print('CHANGEEEEEEEEEEEEEEEE', level)
            print('listNameImgs', listNameImgs)
            print('listNumImgs', listNumImgs)
            print('mdToReduce', mdToReduce)
            print('labelMaxClass', labelMaxClass)

            self._params = {'input': mdToReduce,
                            'numRef': labelMaxClass,
                            #'outputMd': mdToReduce[0:-4]+'_NoAtt.xmd'}
                            'outputMd': mdToReduce[0:-25] + 'test_level00' + str(level) + '_major.xmd'}
            args = ('-i %(input)s --query select "ref==%(numRef)d" -o %(outputMd)s')
            self.runJob("xmipp_metadata_utilities", args % self._params)

        return change, listNumImgs, listNameImgs

            #self.imgsExp = self._getExtraPath('general_test_level00' + str(level) + '_NoAtt.xmd')

            #for i in range(self.numberOfIterations):
            #    self.iterationStep(self.refSet, self.imgsExp, i, False, level, True)
            #    self.refSet = self._getExtraPath('level00' + str(level) + '_classes.xmd')

            #self.attractionSplitStep(level)

            #if self.depth>1:
            #    self.depth-=1
            #    return
            #if self.depth==1:
            #    self.depth=0

            #self.imgsExp = self._getExtraPath('imagesExp.xmd')

            #if level==0:
            #    self.ref = self._getExtraPath('level00' + str(level) + '_classes.xmd')
            #else:
            #    self.generateMetadata(listNameImgs, listNumImgs, level)
            #   self.ref = self._getExtraPath('intermediate_classes.xmd')

            #estas doss lineas no estaban y el siguiente for estaba dentro del while
            #copy(self._getExtraPath('intermediate_classes.xmd'), self._getExtraPath('general_level00' + str(level) + '_classes.xmd'))

            #change, labelMaxClass, labelMinClass, mdToReduce, listNumImgs, listNameImgs = self.checkAttraction(level, False)


        #for i in range(self.numberOfIterations):
        #    self.iterationStep(self.ref, self.exp, i, False, level, False)
        #    self.ref = self._getExtraPath('general_level00' + str(level) + '_classes.xmd')




    def checkAttraction(self, level, flag_split):

        if flag_split:
            mdToCheck = self._getExtraPath('level00'+str(level)+'_classes.xmd')
            mdToReduce = self._getExtraPath('test_level00'+str(level)+'.xmd')
        else:
            mdToCheck = self._getExtraPath('general_level00' + str(level) + '_classes.xmd')
            mdToReduce = self._getExtraPath('general_test_level00' + str(level) + '.xmd')

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

        p = 0.35
        th = (p*total/len(listAuxNum))
        labelMinClass=[]
        print('TH', p * total / len(listAuxNum))
        for i in range(len(listAuxNum)):
            if listAuxNum[i]<th:
                labelMinClass.append(listAuxRef[i])
                listAuxNum[i] = -1

        labelMaxClass = listAuxNum.index(max(listAuxNum))
        listAuxNum[labelMaxClass] = -1

        if len(labelMinClass)>0:
            change = True
        else:
            change = False

        return change, labelMaxClass, labelMinClass, mdToReduce, listAuxNum, listAuxName


    # def checkOutput(self, iter, listNumImgs, count):
    #
    #     if -1 in listNumImgs:
    #         ind = listNumImgs.index(-1)
    #     else:
    #         ind = None
    #     maxValue = 0
    #     listAux = []
    #
    #     outSet = self._getExtraPath('level00'+str(count)+'_iter'+str(iter)+'_classes.xmd')
    #     metadataItem = md.MetaData(outSet)
    #     for item in metadataItem:
    #         numImgs = metadataItem.getValue(md.MDL_CLASS_COUNT, item)
    #         listAux.append(numImgs)
    #
    #     if ind is None:
    #         listNumImgs.append(listAux[0])
    #         listNumImgs.append(listAux[1])
    #         self.levels.append(count)
    #         self.levels.append(count)
    #         self.numRefs.append(0)
    #         self.numRefs.append(1)
    #         #print("listNumImgs ", listNumImgs)
    #     else:
    #         listNumImgs[ind] = listAux[0]
    #         listNumImgs.insert(ind+1, listAux[1])
    #         self.levels[ind]=count
    #         self.levels.insert(ind+1, count)
    #         self.numRefs[ind]=0
    #         self.numRefs.insert(ind+1, 1)
    #         #print("listNumImgs ", listNumImgs)
    #
    #     if len(listNumImgs) is self.numberOfClasses.get():
    #         return
    #
    #     for i in range(len(listNumImgs)):
    #         value = listNumImgs[i]
    #         if value>maxValue:
    #             maxValue = value
    #             maxPos = i
    #     listNumImgs[maxPos] = -1
    #
    #     #print("iter ", iter)
    #     #print("listAux ", listAux)
    #     #print("maxValue ", maxValue)
    #     #print("maxPos ", maxPos)
    #
    #     #print("self.levels ", self.levels)
    #     #print("self.numRefs ", self.numRefs)
    #
    #     # Fifth step: generate a metadata with the most numerous class
    #     testMd = self._getExtraPath('test_level00'+str(self.levels[maxPos])+'.xmd')
    #     self._params = {'input': testMd,
    #                     'numRef': self.numRefs[maxPos],
    #                     'outputMd': self._getExtraPath('test_level00'+str(count)+'_class.xmd')}
    #     args = ('-i %(input)s --query select "ref==%(numRef)d" -o %(outputMd)s')
    #     self.runJob("xmipp_metadata_utilities", args % self._params)



    def checkOutput(self, level):

        print('checkOutput')

        maxValue = 0
        listAuxString = []
        listAuxNum = []
        listAuxRefs = []

        outSet = self._getExtraPath('general_level00'+str(level)+'_classes.xmd')
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

        for i in range(len(listAuxNum)):
            value = listAuxNum[i]
            if value>maxValue:
                maxValue = value
                maxPos = i

        listAuxNum[maxPos] = -1
        bestRef = listAuxRefs[maxPos]
        # name = listAuxString[maxPos]
        # levelImg = name[name.find('level')+5:name.find('level')+8]

        # Fifth step: generate a metadata with the most numerous class
        # testMd = self._getExtraPath('general_test_level'+levelImg+'.xmd')
        # self._params = {'input': testMd,
        #                 'numRef': bestRef,
        #                 'outputMd': self._getExtraPath('test_level00'+str(level)+'_major.xmd')}
        # args = ('-i %(input)s --query select "ref==%(numRef)d" -o %(outputMd)s')
        # self.runJob("xmipp_metadata_utilities", args % self._params)

        mdAll = md.MetaData()
        outSet = self._getExtraPath('general_level00'+str(level)+'_classes.xmd')
        blocks = md.getBlocksInMetaDataFile(outSet)
        for block in blocks:
            if block.startswith('class%06d' % (bestRef)):
                mdClass = md.MetaData(block + "@" + outSet)
                rows = iterRows(mdClass)
                for row in rows:
                    row.addToMd(mdAll)
        mdAll.write(self._getExtraPath('test_level00'+str(level)+'_major.xmd'), MD_APPEND)

        return listAuxNum, listAuxString



    def createOutputStep(self):
        """ Store the setOfParticles object
        as result of the protocol.
        """
        print('createOutputStep')
        #filename = self._getExtraPath('general_level00' + str(self.numberOfClasses.get()-2) + '_classes.xmd')
        filename = self._getExtraPath('last_classes.xmd')
        print(filename)

        imgSet = self.inputParticles.get()
        classes2D = self._createSetOfClasses2D(imgSet)

        readSetOfClasses2D(classes2D, filename)

        self._defineOutputs(outputClasses=classes2D)
        self._defineSourceRelation(self.inputParticles, classes2D)


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

