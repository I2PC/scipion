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


# TODO: finish copying stuff from test_prot_sets2.py and remove this
# comment.

class TestSets(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset_xmipp = DataSet.getDataSet('xmipp_tutorial')
        cls.dataset_mda = DataSet.getDataSet('mda')
        cls.dataset_ribo = DataSet.getDataSet('ribo_movies')

    def testImports(self):
        # Micrographs
        from pyworkflow.em.protocol.protocol_import import ProtImportMicrographs
        p_imp_micros = self.proj.newProtocol(ProtImportMicrographs,
                                            pattern=self.dataset_xmipp.getFile('allMics'),
                                            samplingRate=1.237, voltage=300)
        self.proj.launchProtocol(p_imp_micros, wait=True)
        
        # Volumes
        from pyworkflow.em.protocol.protocol_import import ProtImportVolumes
        p_imp_volumes = self.proj.newProtocol(ProtImportVolumes,
                                         pattern=self.dataset_xmipp.getFile('volumes'),
                                         samplingRate=9.896)
        self.proj.launchProtocol(p_imp_volumes, wait=True)

        # Movies
        from pyworkflow.em.protocol.protocol_import import ProtImportMovies
        p_imp_movies = self.proj.newProtocol(ProtImportMovies,
                                             pattern=self.dataset_ribo.getFile('movies'),
                                             samplingRate=2.37, magnification=59000,
                                             voltage=300, sphericalAberration=2.0)
        self.proj.launchProtocol(p_imp_movies, wait=True)

        # Particles
        from pyworkflow.em.protocol.protocol_import import ProtImportParticles
        p_imp_particles = self.proj.newProtocol(ProtImportParticles,
                                                pattern=self.dataset_mda.getFile('particles'),
                                                samplingRate=3.5)
        self.proj.launchProtocol(p_imp_particles, wait=True)

        # Coordinates
        # Oh, I don't know of any example of coordinates imported :(

    def testSets(self):
        pass


if __name__ == '__main__':
    unittest.main()        
