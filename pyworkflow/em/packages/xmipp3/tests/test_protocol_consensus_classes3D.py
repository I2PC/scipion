#!/usr/bin/env python

# ***************************************************************************
# *
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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

import random

from pyworkflow.tests import BaseTest, setupTestProject, DataSet

from pyworkflow.utils import redStr, greenStr, magentaStr
from pyworkflow.em.data import Class3D

from pyworkflow.em.protocol.protocol_import import ProtImportParticles
from pyworkflow.em.packages.xmipp3.protocol_consensus_classes3D import \
                                                     XmippProtConsensusClasses3D

class TestConsensusClasses3D(BaseTest):

    @classmethod
    def setUpClass(cls):
        """Prepare the data that we will use later on."""

        print "\n", greenStr(" Set Up - Collect data ".center(75, '-'))

        setupTestProject(cls)  # defined in BaseTest, creates cls.proj

        cls.xmippDataTest = DataSet.getDataSet('xmipp_tutorial')


        def _importClasses(label, numClasses, randomSeed, numPart=None):
            """ We import a set of classes with an alterate import particles
                because there is not a import classes protocol
            """
            pImpClasses = cls.proj.newProtocol(ProtImportParticles,
                               filesPath=cls.xmippDataTest.getFile('particles'),
                               samplingRate=3.5)

            pImpClasses.setObjLabel('Import %s' % label)

            # we launch the protocol to obtain the particles
            cls.proj.launchProtocol(pImpClasses, wait=True)

            # fake importing classes handmade
            partSet = pImpClasses.outputParticles

            setOfClasses = pImpClasses._createSetOfClasses3D(partSet)

            numOfPart = partSet.getSize() if numPart is None else numPart
            partIds = list(partSet.getIdSet())
            m = numOfPart/numClasses
            
            # random shuffle with a certain seed to get always the same classes
            random.seed(randomSeed)
            random.shuffle(partIds)
            for clInx in range(0, numClasses):
                currIds = partIds[clInx*m:(clInx+1)*m]

                newClass = Class3D()
                newClass.copyInfo(partSet)
                newClass.setAcquisition(partSet.getAcquisition())
                # newClass.setRepresentative(clRep.getRepresentative())

                setOfClasses.append(newClass)

                enabledClass = setOfClasses[newClass.getObjId()]
                enabledClass.enableAppend()
                for itemId in currIds:
                    item = partSet[itemId]
                    enabledClass.append(item)

                setOfClasses.update(enabledClass)

            # including the outputClasses as protocol output by hand 
            pImpClasses._defineOutputs(outputClasses=setOfClasses)

            # relaunching the protocol to get the new output
            cls.proj.launchProtocol(pImpClasses, wait=True)

            return pImpClasses.outputClasses

        cls.set1 = _importClasses('3 Classes I', 3, 65)
        cls.set2 = _importClasses('3 Classes II', 3, 123)
        cls.set3 = _importClasses('3 Classes III', 3, 256)
        cls.set4 = _importClasses('3 Classes IV (not all part)', 3, 568, 70)
        cls.set5 = _importClasses('5 Classes I', 5, 745)
        cls.set6 = _importClasses('5 Classes II', 5, 1025)


    def checkIntersections(self, setOfIntersections, classId, partIds):
        """ Some assertions for a certain class in the setOfIntersections
        """
        classItem = setOfIntersections[classId]
        self.assertEqual(classItem.getSize(), len(partIds), 
                         "The size of the class %d is wrong" % classId)
        self.assertEqual(partIds,list(classItem.getIdSet()),
                         "The particles in the class %d are wrong." % classId)

    def checkPopulation(self, setOfIntersections, population):
        numOfPart = 0
        for classItem in setOfIntersections:
            numOfPart += classItem.getSize()

        self.assertEqual(numOfPart, population,
                         "The total number of particles is wrong")


    # The tests themselves.
    #
    def testConsensus(self):
        print "\n", greenStr(" Test Consensus ".center(75, '-'))

        # preparing and launching the protocol
        pConsClass = self.proj.newProtocol(XmippProtConsensusClasses3D,
                                           inputMultiClasses=[self.set1, 
                                                              self.set2, 
                                                              self.set3])
        self.proj.launchProtocol(pConsClass, wait=True)
        setOfIntersections = pConsClass.outputClasses

        # some general assertions
        self.assertIsNotNone(setOfIntersections,
                             "There was some problem with the output")
        self.assertEqual(setOfIntersections.getSize(), 27,
                         "The number of the outputClasses is wrong")
        self.checkPopulation(setOfIntersections, 75)

        # some specific assetions 
        self.checkIntersections(setOfIntersections, 1, [53, 54])
        self.checkIntersections(setOfIntersections, 7, [47, 12, 70, 31])
        self.checkIntersections(setOfIntersections, 20,[66, 67, 55])

    def testConsensus2(self):
        print "\n", greenStr(" Test Consensus with different set sizes".center(75, '-'))

        # preparing and launching the protocol
        pConsClass = self.proj.newProtocol(XmippProtConsensusClasses3D,
                                           inputMultiClasses=[self.set1, 
                                                              self.set2, 
                                                              self.set4])
        self.proj.launchProtocol(pConsClass, wait=True)
        setOfIntersections = pConsClass.outputClasses

        # some general assertions
        self.assertIsNotNone(setOfIntersections,
                             "There was some problem with the output")
        self.assertEqual(setOfIntersections.getSize(), 27,
                         "The number of the outputClasses is wrong")
        self.checkPopulation(setOfIntersections, 69)

        # some specific assetions 
        self.checkIntersections(setOfIntersections, 1, [74, 53])
        self.checkIntersections(setOfIntersections, 7, [32, 47, 12, 31])
        self.checkIntersections(setOfIntersections, 20,[48, 67])

    def testConsensus3(self):
        print "\n", greenStr(" Test Consensus with different set sizes 2".center(75, '-'))

        # preparing and launching the protocol
        pConsClass = self.proj.newProtocol(XmippProtConsensusClasses3D,
                                           inputMultiClasses=[self.set4, 
                                                              self.set2, 
                                                              self.set1])
        self.proj.launchProtocol(pConsClass, wait=True)
        setOfIntersections = pConsClass.outputClasses

        # some general assertions
        self.assertIsNotNone(setOfIntersections,
                             "There was some problem with the output")
        self.assertEqual(setOfIntersections.getSize(), 27,
                         "The number of the outputClasses is wrong")
        self.checkPopulation(setOfIntersections, 69)

        # some specific assetions 
        self.checkIntersections(setOfIntersections, 1, [74, 53])
        self.checkIntersections(setOfIntersections, 7, [25, 43, 68, 54, 57])
        self.checkIntersections(setOfIntersections, 20,[46])

    def testConsensus4(self):
        print "\n", greenStr(" Test Consensus with different class sizes".center(75, '-'))

        # preparing and launching the protocol
        pConsClass = self.proj.newProtocol(XmippProtConsensusClasses3D,
                                           inputMultiClasses=[self.set5, 
                                                              self.set2, 
                                                              self.set1])
        self.proj.launchProtocol(pConsClass, wait=True)
        setOfIntersections = pConsClass.outputClasses

        # some general assertions
        self.assertIsNotNone(setOfIntersections,
                             "There was some problem with the output")
        self.assertEqual(setOfIntersections.getSize(), 45,
                         "The number of the outputClasses is wrong")
        self.checkPopulation(setOfIntersections, 75)

        # some specific assetions 
        self.checkIntersections(setOfIntersections, 1, [74, 53])
        self.checkIntersections(setOfIntersections, 7, [25, 57])
        self.checkIntersections(setOfIntersections, 20,[56, 33, 49])

    def testConsensus5(self):
        print "\n", greenStr(" Test Consensus with two sets".center(75, '-'))

        # preparing and launching the protocol
        pConsClass = self.proj.newProtocol(XmippProtConsensusClasses3D,
                                           inputMultiClasses=[self.set5, 
                                                              self.set6])
        self.proj.launchProtocol(pConsClass, wait=True)
        setOfIntersections = pConsClass.outputClasses

        # some general assertions
        self.assertIsNotNone(setOfIntersections,
                             "There was some problem with the output")
        self.assertEqual(setOfIntersections.getSize(), 25,
                         "The number of the outputClasses is wrong")
        self.checkPopulation(setOfIntersections, 75)

        # some specific assetions 
        self.checkIntersections(setOfIntersections, 1, [48, 74, 44, 53])
        self.checkIntersections(setOfIntersections, 7, [51, 68, 14])
        self.checkIntersections(setOfIntersections, 20,[36, 5])

    def testConsensus6(self):
        print "\n", greenStr(" Test Consensus with four sets".center(75, '-'))

        # preparing and launching the protocol
        pConsClass = self.proj.newProtocol(XmippProtConsensusClasses3D,
                                           inputMultiClasses=[self.set1,
                                                              self.set2,
                                                              self.set5, 
                                                              self.set6])
        self.proj.launchProtocol(pConsClass, wait=True)
        setOfIntersections = pConsClass.outputClasses

        # some general assertions
        self.assertIsNotNone(setOfIntersections,
                             "There was some problem with the output")
        self.assertEqual(setOfIntersections.getSize(), 225,
                         "The number of the outputClasses is wrong")
        self.checkPopulation(setOfIntersections, 75)

        # some specific assetions 
        self.checkIntersections(setOfIntersections, 1, [74, 53])
        self.checkIntersections(setOfIntersections, 7, [44])
