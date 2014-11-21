#!/usr/bin/env python

# **************************************************************************
# *
# * Authors: Jordi Burguet Castell (jburguet@cnb.csic.es)
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet

from pyworkflow.em.data import EMObject

from pyworkflow.em.protocol.protocol_import import (
    ProtImportMicrographs, ProtImportVolumes, ProtImportMovies,
    ProtImportParticles, ProtImportCoordinates)

from pyworkflow.em.protocol.protocol_sets import (
    ProtSplitSet, ProtSubSet, ProtUnionSet)



class TestSets(BaseTest):

    """Run different tests related to the set operations."""

    @classmethod
    def setUpClass(cls):
        """Prepare the data that we will use later on."""

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
        p_imp_micros = new(ProtImportMicrographs,
                           filesPath=cls.dataset_xmipp.getFile('allMics'),
                           samplingRate=1.237, voltage=300)
        launch(p_imp_micros, wait=True)
        cls.micros = p_imp_micros.outputMicrographs

        # Volumes
        p_imp_volumes = new(ProtImportVolumes,
                            filesPath=cls.dataset_xmipp.getFile('volumes'),
                            samplingRate=9.896)
        launch(p_imp_volumes, wait=True)
        cls.vols = p_imp_volumes.outputVolumes

        # Movies
        p_imp_movies = new(ProtImportMovies,
                           filesPath=cls.dataset_ribo.getFile('movies'),
                           samplingRate=2.37, magnification=59000,
                           voltage=300, sphericalAberration=2.0)
        launch(p_imp_movies, wait=True)
        cls.movies = p_imp_movies.outputMovies

        # Particles
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

        p_split = self.proj.newProtocol(ProtSplitSet)
        p_split.inputSet.set(em_set)
        p_split.numberOfSets.set(n)
        p_split.randomize.set(randomize)
        self.proj.launchProtocol(p_split, wait=True)
        return p_split

    def outputs(self, p):
        """Iterator over all the elements in the outputs of protocol p."""

        for key, output in p.iterOutputAttributes(EMObject):
            yield output

    #
    # The tests themselves.
    #
    def testSplit(self):
        """Test that the split operation works as expected."""

        def check(set0, n=2, randomize=False):
            "Simple checks on split sets from set0."
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

        def check(set0, n1=2, n2=2):
            "Simple checks on subsets, coming from split sets of set0."
            p_split1 = self.split(set0, n=n1, randomize=True)
            p_split2 = self.split(set0, n=n2, randomize=True)

            setFull = random.choice(list(self.outputs(p_split1)))
            setSub = random.choice(list(self.outputs(p_split2)))
            p_subset = self.proj.newProtocol(ProtSubSet)
            p_subset.inputFullSet.set(setFull)
            p_subset.inputSubSet.set(setSub)
            self.proj.launchProtocol(p_subset, wait=True)

            setFullIds = [x.strId() for x in setFull]
            output = self.outputs(p_subset).next()  # first (and only!) output
            for elem in output:
                self.assertTrue(elem.strId() in setFullIds)

            self.assertTrue(len(setFull) >= len(output))

        # We won't do these first two, there are too few elements.
        #   check(self.micros)
        #   check(self.vols)
        check(self.movies)
        check(self.particles)
        check(self.particles, n1=3, n2=5)

    def testMerge(self):
        """Test that the union operation works as expected."""

        def check(set0):
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



if __name__ == '__main__':
    unittest.main()
