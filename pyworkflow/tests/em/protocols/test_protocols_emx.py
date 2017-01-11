# ***************************************************************************
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
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
# ***************************************************************************/

import os
from itertools import izip

import pyworkflow.tests as tests
from pyworkflow.em.data import SetOfCoordinates, SetOfMicrographs
from pyworkflow.em.protocol import ProtImportMicrographs, ProtImportParticles
from pyworkflow.em.constants import ALIGN_2D
from pyworkflow.em.packages.emxlib import ProtEmxExport


class TestEmxBase(tests.BaseTest):
    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = tests.DataSet.getDataSet('emx')
        cls.dsRelion = tests.DataSet.getDataSet('relion_tutorial')

    def test_emx_export(self):

        prot1 = self.newProtocol(ProtImportParticles,
                                 objLabel='from scipion (to-reconstruct)',
                                 importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                 sqliteFile=self.dsRelion.getFile('import/case2/particles.sqlite'),
                                 magnification=10000,
                                 samplingRate=7.08
                                 )
        self.launchProtocol(prot1)

        protEmx = self.newProtocol(ProtEmxExport)
        protEmx.inputSet.set(prot1.outputParticles)
        self.launchProtocol(protEmx)

    def test_coodinatesTest1(self):
        """ Import an EMX file with just one micrograph
        and a few coordinates.
        """
        protEmxImport   = self.newProtocol(ProtImportParticles,
                                           objLabel='from emx (coordinatesT1)',
                                           importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                           emxFile=self.dataset.getFile('coordinatesT1'),
                                           alignType=3,
                                           voltage=100,
                                           magnification=10000,
                                           samplingRate=2.46)
        self.launchProtocol(protEmxImport)

        # Reference coordinates
        coords = SetOfCoordinates(filename=self.dataset.getFile('coordinatesGoldT1'))
        tests.BaseTest.compareSets(self, protEmxImport.outputCoordinates, coords)

    def test_particleImportDefocus(self):
        """ Import an EMX file with a stack of particles
        that has defocus
        """
        emxFn = self.dataset.getFile('defocusParticleT2')
        protEmxImport = self.newProtocol(ProtImportParticles,
                                         objLabel='emx - import ctf',
                                         importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                         emxFile=emxFn,
                                         alignType=3,
                                         magnification=10000,
                                         samplingRate=2.8
                                         )
        self.launchProtocol(protEmxImport)
        micFn = self.dataset.getFile('micrographsGoldT2')
        mics = SetOfMicrographs(filename = micFn)

        for mic1, mic2 in izip(mics, protEmxImport.outputMicrographs):
            # Remove the absolute path in the micrographs to 
            # really check that the attributes should be equal
            mic1.setFileName(os.path.basename(mic1.getFileName()))
            mic2.setFileName(os.path.basename(mic2.getFileName()))
            self.assertTrue(mic1.equalAttributes(mic2, verbose=True))

    def test_micrographImport(self):
        """ Import an EMX file with micrographs and defocus
        """
        emxFn = self.dataset.getFile('emxMicrographCtf1')
        protEmxImport = self.newProtocol(ProtImportMicrographs,
                                         objLabel='emx - import mics',
                                         importFrom=ProtImportMicrographs.IMPORT_FROM_EMX,
                                         emxFile=emxFn,
                                         magnification=10000,
                                         samplingRate=2.46,
                                         )
        self.launchProtocol(protEmxImport)
        micFn =self.dataset.getFile('emxMicrographCtf1Gold')
        mics  = SetOfMicrographs(filename = micFn)

        for mic1, mic2 in izip(mics, protEmxImport.outputMicrographs):
            # Remove the absolute path in the micrographs to
            # really check that the attributes should be equal
            mic1.setFileName(os.path.basename(mic1.getFileName()))
            mic2.setFileName(os.path.basename(mic2.getFileName()))
            self.assertTrue(mic1.equalAttributes(mic2, verbose=True))


    def test_particlesAlignment(self):
        """ Import an EMX file with particles and alignment information.
        """
        emxFn = self.dataset.getFile('alignment/Test1/images.emx')
        protEmxImport = self.newProtocol(ProtImportParticles,
                                         objLabel='emx: import alignment',
                                         importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                         alignType=0,#ALIGN2D
                                         emxFile=emxFn,
                                         samplingRate=1.,
                                         amplitudeContrast=2.,
                                         sphericalAberration=0.1,
                                         voltage=200)
        self.launchProtocol(protEmxImport)
        self.assertIsNotNone(protEmxImport.outputParticles, "There was a problem with the 'emx: import alignment' outputParticles")
            
    def test_particlesReconstruction(self):
        """ Import an EMX file with particles and alignment information.
        """
        emxFn = self.dataset.getFile('alignment/Test2/stack2D.emx')
        protEmxImport = self.newProtocol(ProtImportParticles,
                                         objLabel='emx: import projection',
                                         importFrom=ProtImportParticles.IMPORT_FROM_EMX,
                                         alignType=2,#PROJ
                                         emxFile=emxFn,
                                         samplingRate=1.,
                                         amplitudeContrast=2.,
                                         sphericalAberration=0.1,
                                         voltage=200)
        self.launchProtocol(protEmxImport)
        self.assertIsNotNone(getattr(protEmxImport, 'outputParticles', None), "There was a problem with the 'emx: import projection' outputParticles")

