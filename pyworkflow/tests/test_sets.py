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


from pyworkflow.tests import BaseTest, setupTestProject, DataSet

from pyworkflow.em.data import EMObject


# TODO: finish copying stuff from test_prot_sets2.py and remove this
# comment.

class TestSets(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset_xmipp = DataSet.getDataSet('xmipp_tutorial')
        cls.dataset_mda = DataSet.getDataSet('mda')
        cls.dataset_ribo = DataSet.getDataSet('ribo_movies')

        #
        # Imports
        #
        new = cls.proj.newProtocol  # short notation
        launch = cls.proj.launchProtocol

        # Micrographs
        from pyworkflow.em.protocol.protocol_import import ProtImportMicrographs
        p_imp_micros = new(ProtImportMicrographs,
                           pattern=cls.dataset_xmipp.getFile('allMics'),
                           samplingRate=1.237, voltage=300)
        launch(p_imp_micros, wait=True)
        cls.micros = p_imp_micros.outputMicrographs

        # Volumes
        from pyworkflow.em.protocol.protocol_import import ProtImportVolumes
        p_imp_volumes = new(ProtImportVolumes,
                            pattern=cls.dataset_xmipp.getFile('volumes'),
                            samplingRate=9.896)
        launch(p_imp_volumes, wait=True)
        cls.vols = p_imp_volumes.outputVolumes

        # Movies
        from pyworkflow.em.protocol.protocol_import import ProtImportMovies
        p_imp_movies = new(ProtImportMovies,
                           pattern=cls.dataset_ribo.getFile('movies'),
                           samplingRate=2.37, magnification=59000,
                           voltage=300, sphericalAberration=2.0)
        launch(p_imp_movies, wait=True)
        cls.movies = p_imp_movies.outputMovies

        # Particles
        from pyworkflow.em.protocol.protocol_import import ProtImportParticles
        p_imp_particles = new(ProtImportParticles,
                              pattern=cls.dataset_mda.getFile('particles'),
                              samplingRate=3.5)
        launch(p_imp_particles, wait=True)
        cls.particles = p_imp_particles.outputParticles

        # Coordinates
        # Oh, I don't know of any example of coordinates imported :(

    def testSplit(self):
        from pyworkflow.em.protocol.protocol_sets import ProtSplitSet

        def split(em_set, n, randomize):
            "Return a run split protocol over input set em_set."
            p_split = self.proj.newProtocol(ProtSplitSet)
            p_split.inputSet.set(em_set)
            p_split.numberOfSets.set(n)
            p_split.randomize.set(randomize)
            self.proj.launchProtocol(p_split, wait=True)
            return p_split

        def outputs(p):
            "Iterate over all the elements in the outputs of protocol p."
            for key, output in p.iterOutputAttributes(EMObject):
                yield output

        def check(set0, n=2, randomize=False):
            split_set = split(set0, n=n, randomize=randomize)
            unsplit_set = [x.strId() for x in set0]
            for em_set in outputs(split_set):
                for elem in em_set:
                    self.assertTrue(elem.strId() in unsplit_set)
            self.assertTrue(sum(len(x) for x in outputs(split_set)) == len(set0))

        check(self.micros)
        check(self.micros, randomize=True)
        check(self.vols)
        check(self.movies)
        check(self.particles)
        check(self.particles, n=4)



if __name__ == '__main__':
    unittest.main()        
