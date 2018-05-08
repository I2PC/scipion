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

from pyworkflow.em.protocol.protocol import EMProtocol
from pyworkflow.protocol.params import MultiPointerParam
from pyworkflow.em import Class3D

class XmippProtConsensusClasses3D(EMProtocol):
    """ Compare several SetOfClasses3D.
        Return the intersection of the input classes.
    """
    _label = 'consensus classes 3D'

    intersectsList = []

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputMultiClasses', MultiPointerParam, important=True,
                      label="Input Classes", pointerClass='SetOfClasses3D',
                      help='Select several sets of classes where '
                           'to evaluate the intersections.')

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        """ Inserting one step for each intersections analisis
        """

        self._insertFunctionStep('compareFirstStep', 
                                 self.inputMultiClasses[0].get().getObjId(),
                                 self.inputMultiClasses[1].get().getObjId())

        if len(self.inputMultiClasses)>2:
            for i in range(2,len(self.inputMultiClasses)):
                self._insertFunctionStep('compareOthersStep', i,
                                     self.inputMultiClasses[i].get().getObjId())

        self._insertFunctionStep('createOutputStep')

    def compareFirstStep(self, objId1, objId2):
        """ We found the intersections for the two firsts sets of classes
        """
        set1Id = 0
        set2Id = 1
        set1 = self.inputMultiClasses[set1Id].get()
        set2 = self.inputMultiClasses[set2Id].get()

        print('Computing intersections between classes form set %s and set %s:'
               % (set1.getNameId(), set2.getNameId()))

        newList = []
        for cls1 in set1:
            cls1Id = cls1.getObjId()
            ids1 = cls1.getIdSet()

            for cls2 in set2:
                cls2Id = cls2.getObjId()
                ids2 = cls2.getIdSet()

                interTuple = self.intersectClasses(set1Id, cls1Id, ids1,
                                                   set2Id, cls2Id, ids2)

                newList.append(interTuple)

        self.intersectsList = newList

    def compareOthersStep(self, set1Id, objId):
        """ We found the intersections for the two firsts sets of classes
        """
        set1 = self.inputMultiClasses[set1Id].get()

        print('Computing intersections between classes form set %s and '
              'the previous ones:' % (set1.getNameId()))

        newList = []
        currDB = self.intersectsList
        for cls1 in set1:
            cls1Id = cls1.getObjId()
            ids1 = cls1.getIdSet()

            for currTuple in currDB:
                ids2 = currTuple[1]
                set2Id = currTuple[2]
                cls2Id = currTuple[3]
                clSize = currTuple[4]

                interTuple = self.intersectClasses(set1Id, cls1Id, ids1,
                                                   set2Id, cls2Id, ids2, clSize)

                newList.append(interTuple)
                
        self.intersectsList = newList

    def createOutputStep(self):

        # self.intersectsList.sort(key=lambda e: e[0], reverse=True)

        # print("   ---   S O R T E D:   ---")
        # numberOfPart = 0
        # for classItem in self.intersectsList:
        #     printTuple = (classItem[0], classItem[2])
        #     print(printTuple)
        #     numberOfPart += classItem[0]

        # print('Number of intersections: %d' % len(self.intersectsList))
        # print('Total of particles: %d' % numberOfPart)

        inputParticles = self.inputMultiClasses[0].get().getImages()
        outputClasses = self._createSetOfClasses3D(inputParticles)

        for classItem in self.intersectsList:
            numOfPart = classItem[0]
            partIds = classItem[1]
            setRepId = classItem[2]
            clsRepId = classItem[3]

            setRep = self.inputMultiClasses[setRepId].get()
            clRep = setRep[clsRepId]

            newClass = Class3D()
            newClass.copyInfo(clRep)
            newClass.setAcquisition(clRep.getAcquisition())
            newClass.setRepresentative(clRep.getRepresentative())

            outputClasses.append(newClass)

            enabledClass = outputClasses[newClass.getObjId()]
            enabledClass.enableAppend()
            for itemId in partIds:
                enabledClass.append(inputParticles[itemId])

            outputClasses.update(enabledClass)

        self._defineOutputs(outputClasses=outputClasses)
        for item in self.inputMultiClasses:
            self._defineSourceRelation(item, outputClasses)


    # --------------------------- INFO functions -------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = [] if len(self.inputMultiClasses)>1 else \
                 ["More than one Input Classes is needed to compute the consensus."]
        return errors


    # --------------------------- UTILS functions ------------------------------
    def intersectClasses(self, setId1, clId1, ids1,
                               setId2, clId2, ids2, clsSize2=None):
        size1 = len(ids1)
        size2 = len(ids2) if clsSize2 is None else clsSize2

        inter = ids1.intersection(ids2)

        if size1 < size2:
            setId = setId1
            clsId = clId1
            clsSize = size1
        else:
            setId = setId2
            clsId = clId2
            clsSize = size2

        # print(" ")
        # print(" - Intersection of cl%d of set%d (%d part.) and "
        #                          "cl%d of set%d (%d part.):"
        #        % (clId1, setId1, len(ids1), clId2, setId2, len(ids2)))
        # print("    Size1=%d < Size2=%d = %s" 
        #        % (size1, size2, size1<size2))
        # print("      -> from set %d calss %d, with %d part. in the intersection." 
        #        % (setId, clsId, len(inter)))
        # print(" -  -  -  -  -  -  -  -  -  -")

        return (len(inter), inter, setId, clsId, clsSize)
