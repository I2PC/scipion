#!/usr/bin/env python

# ***************************************************************************
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Jordi Burguet Castell (jburguet@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import random
import unittest
from itertools import izip

from pyworkflow.tests import BaseTest, setupTestProject, DataSet

from pyworkflow.utils import redStr, greenStr, magentaStr
from pyworkflow.em.data import EMObject

from pyworkflow.em.protocol.protocol_import import (
    ProtImportMicrographs, ProtImportVolumes, ProtImportMovies,
    ProtImportParticles, ProtImportCoordinates)

from pyworkflow.em.protocol.protocol_sets import (
    ProtSplitSet, ProtSubSet, ProtUnionSet)

# Used by Roberto's test, where he creates the particles "by hand"
from pyworkflow.em.data import Particle, SetOfParticles, Acquisition
from pyworkflow.utils.utils import prettyDict


class TestSets(BaseTest):

    """Run different tests related to the set operations."""

    @classmethod
    def setUpClass(cls):
        """Prepare the data that we will use later on."""

        print "\n", greenStr(" Set Up - Collect data ".center(75, '-'))

        setupTestProject(cls)  # defined in BaseTest, creates cls.proj

        cls.dataset_xmipp = DataSet.getDataSet('xmipp_tutorial')
        cls.dataset_mda = DataSet.getDataSet('mda')
        cls.dataset_ribo = DataSet.getDataSet('ribo_movies')

        #
        # Imports
        #
        new = cls.proj.newProtocol  # short notation
        launch = cls.proj.launchProtocol
        # Micrographs
        print magentaStr("\n==> Importing data - micrographs")
        p_imp_micros = new(ProtImportMicrographs,
                           filesPath=cls.dataset_xmipp.getFile('allMics'),
                           samplingRate=1.237, voltage=300)
        launch(p_imp_micros, wait=True)
        cls.micros = p_imp_micros.outputMicrographs

        # Volumes
        print magentaStr("\n==> Importing data - volumes")
        p_imp_volumes = new(ProtImportVolumes,
                            filesPath=cls.dataset_xmipp.getFile('volumes'),
                            samplingRate=9.896)
        launch(p_imp_volumes, wait=True)
        cls.vols = p_imp_volumes.outputVolumes

        # Movies
        print magentaStr("\n==> Importing data - movies")
        p_imp_movies = new(ProtImportMovies,
                           filesPath=cls.dataset_ribo.getFile('movies'),
                           samplingRate=2.37, magnification=59000,
                           voltage=300, sphericalAberration=2.0)
        launch(p_imp_movies, wait=True)
        cls.movies = p_imp_movies.outputMovies

        # Particles
        print magentaStr("\n==> Importing data - particles")
        p_imp_particles = new(ProtImportParticles,
                              filesPath=cls.dataset_mda.getFile('particles'),
                              samplingRate=3.5)
        launch(p_imp_particles, wait=True)
        cls.particles = p_imp_particles.outputParticles

        # Coordinates
        # Oh, I don't know of any example of coordinates imported :(

    #
    # Helper functions
    #
    def split(self, em_set, n, randomize):
        """Return a run split protocol over input set em_set."""

        p_split = self.proj.newProtocol(ProtSplitSet, numberOfSets=n)
        p_split.inputSet.set(em_set)
        p_split.randomize.set(randomize)
        self.proj.launchProtocol(p_split, wait=True)
        return p_split

    def outputs(self, p):
        """Iterator over all the elements in the outputs of protocol p."""
        for key, output in p.iterOutputEM():
            yield output
    #
    # The tests themselves.
    #
    def testSplit(self):
        """Test that the split operation works as expected."""

        print "\n", greenStr(" Test Split ".center(75, '-'))

        def check(set0, n=2, randomize=False):
            "Simple checks on split sets from set0."
            print magentaStr("\n==> Check split of %s" % type(set0).__name__)
            unsplit_set = [x.strId() for x in set0]
            p_split = self.split(set0, n=n, randomize=randomize)
            # Are all output elements of the protocol in the original set?
            for em_set in self.outputs(p_split):
                for elem in em_set:
                    self.assertTrue(elem.strId() in unsplit_set)
            # Number of elements of all splitted sets equal to original number?
            self.assertEqual(sum(len(x) for x in self.outputs(p_split)),
                             len(set0))

        check(self.micros)
        check(self.micros, randomize=True)
        check(self.vols)
        check(self.movies)
        check(self.particles)
        check(self.particles, n=4)

    def testSubset(self):
        """Test that the subset operation works as expected."""

        print "\n", greenStr(" Test Subset ".center(75, '-'))

        def check(set0, n1=2, n2=2):
            "Simple checks on subsets, coming from split sets of set0."
            print magentaStr("\n==> Check subset of %s" % type(set0).__name__)
            p_split1 = self.split(set0, n=n1, randomize=True)
            p_split2 = self.split(set0, n=n2, randomize=True)

            setFull = random.choice(list(self.outputs(p_split1)))
            setSub = random.choice(list(self.outputs(p_split2)))
            
            label = '%s - %s,%s ' % (set0.getClassName(), n1, n2)
            # Launch intersection subset
            p_subset = self.newProtocol(ProtSubSet)
            p_subset.setObjLabel(label + 'intersection')
            p_subset.inputFullSet.set(setFull)
            p_subset.inputSubSet.set(setSub)
            self.launchProtocol(p_subset)
            
            # Launch difference subset
            p_subset_diff = self.proj.copyProtocol(p_subset)
            p_subset_diff.setOperation.set(p_subset_diff.SET_DIFFERENCE)
            p_subset_diff.setObjLabel(label + 'difference')
            self.launchProtocol(p_subset_diff)

            setFullIds = setFull.getIdSet()
            setSubIds = setSub.getIdSet()
            n = len(setFull)
            
            # Check intersection
            outputs = [o for o in self.outputs(p_subset)]
            n1 = 0
            if outputs:
                output = outputs[0]
                n1 = len(output)
                for elem in output:
                    self.assertTrue(elem.getObjId() in setFullIds)
                    self.assertTrue(elem.getObjId() in setSubIds,
                                    'object id %s not in set: %s' % (elem.getObjId(), setSubIds))
                
            # Check difference
            outputs = [o for o in self.outputs(p_subset_diff)]
            n2 = 0
            if outputs:
                output_diff = outputs[0]
                n2 = len(output_diff)
                for elem in output_diff:
                    self.assertTrue(elem.getObjId() in setFullIds)
                    self.assertTrue(elem.getObjId() not in setSubIds)
                
                
            self.assertTrue(n >= n1)
            self.assertTrue(n >= n2)            
            self.assertEqual(n, n1+n2)

        # We won't do these first two, there are too few elements.
        #   check(self.micros)
        #   check(self.vols)
        check(self.movies)
        check(self.particles)
        check(self.particles, n1=3, n2=5)

    def testMerge(self):
        """Test that the union operation works as expected."""

        print "\n", greenStr(" Test Merge ".center(75, '-'))

        def check(set0):
            "Simple checks on merge, coming from many split sets of set0."
            print magentaStr("\n==> Check merge of %s" % type(set0).__name__)
            p_union = self.proj.newProtocol(ProtUnionSet)

            setsIds = []
            for i in range(random.randint(1, 5)):
                n = random.randint(1, len(set0) // 2)
                p_split = self.split(set0, n=n, randomize=True)
                setRandom = random.choice(list(self.outputs(p_split)))
                setsIds.append([x.strId() for x in setRandom])
                p_union.inputSets.append(setRandom)
            self.proj.launchProtocol(p_union, wait=True)

            output = self.outputs(p_union).next()  # first (and only!) output
            self.assertEqual(len(output), sum(len(x) for x in setsIds))
            # We might be able to do more interesting tests, using the
            # collected setsIds.

        check(self.micros)
        check(self.vols)
        check(self.movies)
        check(self.particles)

    def testMergeDifferentAttrs(self):
        """ Test merge from subsets with different attritubes. """
        partFn = self.dataset_mda.getFile('particles/xmipp_particles.xmd')
        protImport = self.newProtocol(ProtImportParticles,
                                      objLabel='import from xmd',
                                      importFrom=ProtImportParticles.IMPORT_FROM_XMIPP3,
                                      mdFile=partFn,
                                      samplingRate=3.5)

        self.launchProtocol(protImport)
        protJoin = self.newProtocol(ProtUnionSet,
                                    objLabel='join different',
                                    ignoreExtraAttributes=True)

        protJoin.inputSets.append(self.particles)
        protJoin.inputSets.append(protImport.outputParticles)

        self.launchProtocol(protJoin)


    def testOrderBy(self):
        """ create set of particles and orderby a given attribute
        """
        # This function was written by Roberto. It does things
        # differently, so let's keep it for reference.

        #create set of particles

        inFileNameMetadata = self.proj.getTmpPath('particlesOrderBy.sqlite')
        inFileNameData = self.proj.getTmpPath('particlesOrderBy.stk')

        imgSet = SetOfParticles(filename=inFileNameMetadata)
        imgSet.setSamplingRate(1.5)
        acq = Acquisition()
        acq.setAmplitudeContrast(0.1)
        acq.setMagnification(10000)
        acq.setVoltage(200)
        acq.setSphericalAberration(2.0)
        
        imgSet.setAcquisition(acq)
        img = Particle()

        for i in range(1, 10):
            img.setLocation(i, inFileNameData)
            img.setMicId(i%3)
            img.setClassId(i%5)
            imgSet.append(img)
            img.cleanObjId()

        imgSet.write()
        #now import the dataset
        prot1 = self.newProtocol(ProtImportParticles,
                                 importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                 sqliteFile=inFileNameMetadata,
                                 magnification=10000,
                                 samplingRate=1.5
                                 )
        prot1.setObjLabel('from sqlite (test-sets)')
        self.launchProtocol(prot1)

        if prot1.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % inFileNameMetadata)
        
        protSplitSet   = self.newProtocol(ProtSplitSet,
                                          inputSet=prot1.outputParticles,
                                          numberOfSets=2,
                                          randomize=True)
        self.launchProtocol(protSplitSet)

        inputSets = [protSplitSet.outputParticles01,protSplitSet.outputParticles02]
        outputSet = SetOfParticles(filename=self.proj.getTmpPath('gold.sqlite'))
        for itemSet in inputSets:
            for obj in itemSet:
                outputSet.append(obj)

        for item1, item2 in izip(imgSet, outputSet):
            if not item1.equalAttributes(item2):
                print "Items differ:"
                prettyDict(item1.getObjDict())
                prettyDict(item2.getObjDict())
            self.assertTrue(item1.equalAttributes(item2),  )



if __name__ == '__main__':
    unittest.main()
