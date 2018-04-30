# **************************************************************************
# *
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es)
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

from pyworkflow.object import Float, Object, String, Float, Integer
from pyworkflow.utils.path import cleanPath
from pyworkflow.protocol.params import MultiPointerParam

from pyworkflow.em import EMSet, Class3D
from pyworkflow.em.protocol.protocol import EMProtocol

# For the viewer part
from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import DataView


class ProtClassesConsensus4(EMProtocol):
    """ Compare several SetOfClasses
    """
    _label = 'consensus 3Dclasses (new 2)'

    interDB = []

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputMultiClasses', MultiPointerParam, important=True,
                      label="Input Classes", pointerClass='SetOfClasses',
                      help='Select several sets of classes where '
                           'to evaluate the consensus.')

        form.addParallelSection(threads=4, mpi=0)

#--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Inserting one step for each couple of setOfClasses comparision
        """

        self._insertFunctionStep('compareFirstStep', 0, 1)

        if len(self.inputMultiClasses)>2:
            for i in range(2,len(self.inputMultiClasses)):
                self._insertFunctionStep('compareOthersStep', i, self.interDB)


        self._insertFunctionStep('createOutputStep')


    def compareFirstStep(self, idxSetCls1, idxSetCls2):

        print('idxSetCls1 = %d'%(idxSetCls1))
        print('idxSetCls2 = %d'%(idxSetCls2))
        
        set1 = self.inputMultiClasses[idxSetCls1].get()
        set2 = self.inputMultiClasses[idxSetCls2].get()

        for cls1 in set1:
            ids1 = cls1.getIdSet()

            for cls2 in set2:
                ids2 = cls2.getIdSet()

                inter = ids1.intersection(ids2)

                if len(ids1) < len(ids2):
                    repCls = cls1
                else:
                    repCls = cls2

                interTuple = (len(inter), inter, repCls)
                self.interDB.append(interTuple)


    def compareOthersStep(self, idxSetCls1, currDB):

        set1 = self.inputMultiClasses[idxSetCls1].get()

        print('idxSetCls1 = %d'%(idxSetCls1))

        newDB = []

        for cls1 in set1:
            ids1 = cls1.getIdSet()

            for cls2 in currDB:
                ids2 = cls2[1]

                inter = ids1.intersection(cls2[1])

                if len(ids1) < len(ids2):
                    repCls = cls1
                else:
                    repCls = cls2[2]

                interTuple = (len(inter), inter, repCls)
                newDB.append(interTuple)
                
        self.interDB = newDB




    def createOutputStep(self):

        # numberOfPart = 0
        # for iii in range(0,len(self.interDB)):
        #
        #     printTuple = (self.interDB[iii][0], self.interDB[iii][2])
        #     print(printTuple)
        #     numberOfPart += self.interDB[iii][0]

        # self.interDB.sort(key=lambda e: e[0], reverse=True)

        # print("   ---   S O R T E D:   ---")
        # numberOfPart = 0
        # for iii in range(0, len(self.interDB)):

        #     printTuple = (self.interDB[iii][0], self.interDB[iii][2])
        #     print(printTuple)
        #     numberOfPart += self.interDB[iii][0]

        # print('Number of intersections: %d' % len(self.interDB))
        # print('Total of particles: %d' % numberOfPart)

        # outputFn = self._getPath('consensus.sqlite')
        # cleanPath(outputFn)
        inputParticles = self.inputMultiClasses[0].get().getImages()
        outputClasses = self._createSetOfClasses3D(inputParticles)
        # classesIds = []

        for clInx in range(0, len(self.interDB)):
            numOfPart = self.interDB[clInx][0]
            partIds = self.interDB[clInx][1]
            clasRep = self.interDB[clInx][2]

            newClass = Class3D()
            newClass.copyInfo(clasRep)
            newClass.setAcquisition(clasRep.getAcquisition())
            newClass.setRepresentative(clasRep.getRepresentative())

            outputClasses.append(newClass)

            enabledClass = outputClasses[newClass.getObjId()]
            # classesIds.append(newClass.getObjId())
            enabledClass.enableAppend()
            for itemId in partIds:
                enabledClass.append(inputParticles[itemId])

            outputClasses.update(enabledClass)

        self._defineOutputs(outputConsensus=outputClasses)
        for i in range(0,len(self.inputMultiClasses)):
            self._defineSourceRelation(self.inputMultiClasses[i], outputClasses)


        # print("   ---   O U T P U T:   ---")
        # numberOfPart = 0
        # print(outputClasses)
        # for iii in classesIds:
        #     classe = outputClasses[iii]
        #     printTuple = (classe.getSize(), classe)
        #     print(printTuple)
        #     numberOfPart += classe.getSize()

        # print('Number of intersections: %d' % len(outputClasses))
        # print('Total of particles: %d' % numberOfPart)


        
        #
        # for clInx in range(0, len(self.interDB)):
        #     partIds = self.interDB[clInx][1]
        #     clasRep = self.interDB[clInx][2]




        #     newClass = Class2D(objId=repId)
        #     newClass.setAlignment2D()
        #     newClass.copyInfo(inputSet)
        #     newClass.setAcquisition(inputSet.getAcquisition())
        #     newClass.setRepresentative(rep)
        #     outputClasses.append(newClass)
        #
        # i = 1
        # mdBlocks = md.getBlocksInMetaDataFile(self._getExtraPath(
        #     'final_classes.xmd'))
        # for block in mdBlocks:
        #     if block.startswith('class00'):
        #         mdClass = md.MetaData(block + "@" + self._getExtraPath(
        #             'final_classes.xmd'))
        #         imgClassId = i
        #         newClass = outputClasses[imgClassId]
        #         newClass.enableAppend()
        #         for row in md.iterRows(mdClass):
        #             part = rowToParticle(row)
        #             newClass.append(part)
        #         i += 1
        #         newClass.setAlignment2D()
        #         outputClasses.update(newClass)
        #






    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        # numPart = self.inputMultiClasses[0].get().getSize()
        # for inClassId in range(0,len(self.inputMultiClasses)):
        #     if self.inputMultiClasses[inClassId].get().getSize() == numPart:
        #         errors.append('All classes should come from the same original '
        #                       'set of particles. Class %d has different size' 
        #                       %inClassId)
        return errors
    
    
# class ViewerClassesConsensus(Viewer):
#     _environments = [DESKTOP_TKINTER, WEB_DJANGO]
#     _targets = [ProtClassesConsensus4]
    
#     def _visualize(self, obj, **kwargs):
#         labels = 'class1.id class1._representative._filename class2.id class2._representative._filename jaccard intersection union'
#         return [DataView(obj.outputConsensus.getFileName(), 
#                          viewParams={'order': labels, 'mode': 'metadata',
#                                      'visible': labels,
#                                      'render': 'class1._representative._filename class2._representative._filename'
#                                      })
#                 ]
    
#     def visualize(self, obj, **kwargs):
#         self._visualize(obj, **kwargs)[0].show()
