# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import sys, unittest

from pyworkflow.em import ProtImportParticles, ProtImportVolumes
from pyworkflow.tests import *
from pyworkflow.em.packages.localized_recons import *


# Some utility functions to import micrographs that are used
# in several tests.
class TestLocalizedReconsBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='relion_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.particlesFn = cls.dataset.getFile('import/refine3d/extra/relion_it001_data.star')
        cls.vol = cls.dataset.getFile('volume')
    
    @classmethod
    def runImportParticles(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportParticles, 
                                     objLabel='from relion (auto-refine 3d)',
                                     importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                     starFile=pattern,
                                     magnification=65000,
                                     samplingRate=samplingRate,
                                     haveDataBeenPhaseFlipped=True
                                      )
        cls.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
        return protImport
    
    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportVolumes, 
                                     filesPath=pattern, samplingRate=samplingRate)
        cls.launchProtocol(protImport)
        return protImport


class TestLocalizedRecons(TestLocalizedReconsBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestLocalizedReconsBase.setData('relion_tutorial')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 7.08)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 7.08)
      
    def testProtLocalizedReconstruction(self):
        print "Run ProtLocalizedReconstruction"
        localizedRecons = self.newProtocol(ProtLocalizedRecons,
                                        splitParticles=True,
                                        boxSize=10,
                                        defineVector=1)
        localizedRecons.inputParticles.set(self.protImport.outputParticles)
        localizedRecons.inputVolume.set(self.protImportVol.outputVolume)
        self.launchProtocol(localizedRecons)
        self.assertIsNotNone(localizedRecons.outputParticles,"There was a "
                             "problem with localized reconstruction protocol")
        
