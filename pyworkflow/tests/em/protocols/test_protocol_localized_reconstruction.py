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
from pyworkflow.em.packages.opic import *


# Some utility functions to import micrographs that are used
# in several tests.
class TestLocalizedReconsBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_programs'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.particlesFn = cls.dataset.getFile('input/Protocol_Projection_Match'
                                              'ing/ProjMatch/goldStandard/Iter'
                                              '_01/current_angles.xmd')
        cls.chimeraFile = cls.dataset.getFile('input/Protocol_Projection_Match'
                                              'ing/ico.cmm')
    
    @classmethod
    def runImportParticles(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportParticles, 
                                     objLabel='from Xmipp ProjMatch',
                                     importFrom=ProtImportParticles.IMPORT_FROM_XMIPP3,
                                     mdFile=pattern,
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
        TestLocalizedReconsBase.setData('xmipp_programs')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 1)
    #         cls.protImportVol = cls.runImportVolumes(cls.vol, 1)

    def _runSubparticles(self, checkSize, **kwargs):
        label = 'localized subpartices ('
        for t in kwargs.iteritems():
            label += '%s=%s' % t
        label += ')'

        prot = self.newProtocol(ProtLocalizedRecons,
                                 objLabel=label,
                                 symmetryGroup="I3",
                                 defineVector=0,
                                 unique=5, **kwargs)
        prot.vectorFile.set(self.chimeraFile)
        prot.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputCoordinates,
                             "There was a problem with localized subparticles protocol")
        self.assertEqual(checkSize, prot.outputCoordinates.getSize())

        return prot

    def testProtLocalizedReconstruction(self):
        print "Run ProtLocalized Reconstruction"

        localSubparticles = self._runSubparticles(120)

        self._runSubparticles(94, mindist=10)
        self._runSubparticles(47, side=25)
        self._runSubparticles(20, top=50)

        localExtraction = self.newProtocol(ProtLocalizedExtraction, boxSize=14)
        localExtraction.inputParticles.set(self.protImport.outputParticles)
        localExtraction.inputCoordinates.set(localSubparticles.outputCoordinates)
        self.launchProtocol(localExtraction)
        self.assertIsNotNone(localExtraction.outputParticles,
                             "There was a problem with localized extraction protocol")
        
        
